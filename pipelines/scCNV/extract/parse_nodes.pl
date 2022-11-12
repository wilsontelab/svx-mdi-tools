use strict;
use warnings;

# constants
use constant {
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNAME_INDEX => 6,
    #-------------
    _IS_PAIRED => 1, # SAM FLAG bits
    _PROPER => 2,
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
    #-------------
    _SHIFT_REVERSE => 4, # how far to shift those bits to yield binary values
    #-------------
    MOL_ID => 0,    # molecule-level data, carried here in QNAME
    UMI1 => 1,      # appended to QNAME by genomex-mdi-tools align / svCapture make_consensus.pl
    UMI2 => 2,      #     will be added by extract_nodes if user-bam was prepared elsewhere
    IS_MERGED => 3,
    IS_DUPLEX => 4, # appended to QNAME by svCapture make_consensus.pl (but not genomex-mdi-tools align)
    STRAND_COUNT1 => 5,
    STRAND_COUNT2 => 6,
    MOL_CLASS => 7, # values added by extract_nodes regardless of pipeline or bam source
    MOL_STRAND => 8,
    IS_OUTER_CLIP1 => 9,
    IS_OUTER_CLIP2 => 10,
    TARGET_CLASS => 11, # values added (or initialized) for svCapture
    SHARED_PROPER => 12, 
    #-------------
    _CLIP => 0, # outData and innData fields
    _POS  => 1,
    _SIDE => 2,
    _SEQ  => 3,
    _NODE => 4,
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    #-------------
    LEFT  => 0, # clip recovery direction for code readability
    RIGHT => 1,   
    #-------------
    LEFTWARD  => "L", # orientation of aligned source molecule relative to a mapped endpoint
    RIGHTWARD => "R",
    #-------------
    IS_PROPER => 'P', # proper and anomalous molecule classes
    IS_SV     => 'V',
    #-------------
    GAP        => 0, # SV evidence type codes, i.e. node classes
    SPLIT      => 1,
    OUTER_CLIP => 2,
    #-------------
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "I",
    DUPLICATION   => "D",
    DELETION      => "L",
    UNKNOWN       => "?",
    PROPER        => 'P', 
    #-------------
    _NULL => '*',   
    #-------------
    MAX_FAMILY_COUNT => 100,
    #-------------
    BIN_SIZE => 20000 # fixed to match 10x scCNV
};

# operating parameters
use vars qw($leftClip_ $rightClip_);
my @clips = ($leftClip_, $rightClip_);
my @sides = (RIGHTWARD,  LEFTWARD); # e.g., at a left end clip, the read continues righward from that point

# working variables
use vars qw(@alns @jxns $molKey
            $maxInsertSize %insertSizes
            $MAX_TLEN $READ_LEN);  
my $gapDelLimit =  2 * $MAX_TLEN - 2 * $READ_LEN; # these are conservative, i.e., call more proper than aligner typically
my $gapDupLimit = -2 * $READ_LEN;

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
# source molecules that are considered proper (merged or unmerged)
sub commitProperMolecule { # aln1/umi1 are the left/proximal end of the source molecule as provided to the aligner
    my ($aln1, $aln2) = @_;
    setEndpointData(\my @outData1, $aln1, LEFT,  RIGHT); # LEFT = process the left clip if mapped on the forward strand
    setEndpointData(\my @outData2, $aln2, RIGHT, LEFT);  
    setOuterEndpoints($aln1, $aln2, \@outData1, \@outData2);  
    # # prepare to write a cross-tabulated file of proper molecule insert sizes
    # my $insertSize = abs($outData2[_POS] - $outData1[_POS]) + 1 + $outData1[_CLIP] + $outData2[_CLIP];
    # $insertSize <= $maxInsertSize and $insertSizes{$mol[TARGET_CLASS]}[$insertSize]++;
}
#---------------------------------------------------------------------------------------------------
# a merged read-pair (or single read) aligned in at least two segments, i.e., that crossed an SV junction
sub parseMergedSplit {
    # shortest clip1 will match umi1, longest will match umi2
    my ($outData1s, @alnIs) = sortReadAlignments(\@alns, LEFT, RIGHT);
    setEndpointData(\my @outData2, $alns[$alnIs[$#alnIs]], RIGHT, LEFT);
    setOuterEndpoints(@alns[$alnIs[0], $alnIs[$#alnIs]], $$outData1s[$alnIs[0]], \@outData2);     
    setSequencedJunctions(@alns[@alnIs]); 
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read pair with no split alignments (i.e., exactly one aln each read)
# but anomalous due to an unaligned junction (possibly at an inner clip)
sub parseUnmergedHiddenJunction {
    my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
    # override aligner proper flag sometimes (aligner fails to call some things proper, e.g., failed-merge overruns)
    setEndpointData(\my @innData1, $alns[$read1], RIGHT, LEFT);
    setEndpointData(\my @innData2, $alns[$read2], LEFT,  RIGHT);
    my $jxnType = getJxnType($alns[$read1], $alns[$read2], \@innData1, \@innData2, GAP);
    $jxnType eq PROPER and return commitProperMolecule($alns[$read1], $alns[$read2], RIGHT, LEFT);    
    # commit SV if truly not proper    
    setEndpointData(\my @outData1, $alns[$read1], LEFT,  RIGHT);
    setEndpointData(\my @outData2, $alns[$read2], RIGHT, LEFT);
    setOuterEndpoints($alns[$read1], $alns[$read2], \@outData1, \@outData2); # actions based on inner portion of molecule    
    setJunction(GAP, $alns[$read1], $alns[$read2]);      
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read-pair that sequenced an SV junction in >=1 of its reads (the most complex alignment patterns land here)
sub parseUnmergedSplit {
    # determine which reads were split
    my (@alnsByRead, @alnIs, @outDatas);
    foreach my $aln(@alns){
        if(!($$aln[FLAG] & _UNMAPPED)){  # unmapped v. rare since remapping; discard them ***********
            my $read = $$aln[FLAG] & _FIRST_IN_PAIR ? READ1: READ2;
            push @{$alnsByRead[$read]}, $aln;
        }
    }     
    # catch complex pairs where one unmerged read didn't get to us
    ($alnsByRead[READ1] and $alnsByRead[READ2]) or return parseMergedSplit();
    # catch prior merge failures at over-running read pairs
    if(@{$alnsByRead[READ1]} == @{$alnsByRead[READ2]}){
        my $key1 = join(":", sort map { join(":", @$_[RNAME_INDEX,POS]) } @{$alnsByRead[READ1]});
        my $key2 = join(":", sort map { join(":", @$_[RNAME_INDEX,POS]) } @{$alnsByRead[READ2]});
        if($key1 eq $key2){ # same set of split alignments for each read in an unmerged pair
            @alns = @{$alnsByRead[READ1]}; # eliminate one read, it is identical to the other
            return parseMergedSplit();
        }  
    }
    # determine the order of the alignments along the sequenced molecule in the split reads
    foreach my $read(READ1, READ2){
        my ($fwdSide, $revSide) = ($read == READ1) ? (LEFT, RIGHT) : (RIGHT, LEFT);
        if(@{$alnsByRead[$read]} == 1){ # an unsplit read paired with a split read
            @{$alnIs[$read]} = (0); # the index within @alnsByRead, not @alns
            setEndpointData(\my @outData, $alnsByRead[$read][0], $fwdSide, $revSide);
            $outDatas[$read] = [\@outData];
        } else { # one of either one or two split reads
            ($outDatas[$read], @{$alnIs[$read]}) = sortReadAlignments($alnsByRead[$read], $fwdSide, $revSide); 
        }  
    }
    # actions based on outer portion of total source molecule
    my $outI1 = $alnIs[READ1][0];
    my $outI2 = $alnIs[READ2][0]; 
    setOuterEndpoints(
        $alnsByRead[READ1][$outI1], $alnsByRead[READ2][$outI2],
          $outDatas[READ1][$outI1],   $outDatas[READ2][$outI2]
    );    
    # actions based on inner portion of molecule
    # all junctions in all split reads
    foreach my $read(READ1, READ2){
        @{$alnsByRead[$read]} > 1 and
            setSequencedJunctions(@{$alnsByRead[$read]}[@{$alnIs[$read]}]);
    }  
    # possible additional SV junction in the gap
    my $innI1 = $alnIs[READ1][$#{$alnIs[READ1]}];
    my $innI2 = $alnIs[READ2][$#{$alnIs[READ2]}]; 
    setJunction(GAP, $alnsByRead[READ1][$innI1], $alnsByRead[READ2][$innI2]); 
} 
#===================================================================================================

#===================================================================================================
# handle SV junction output
#---------------------------------------------------------------------------------------------------
# a read definitively crosses a junction via one or more split reads
sub setSequencedJunctions{ 
    my (@alns) = @_; # two or more alignments for a single contiguous (perhaps merged) read
    
    # reorder aln2 splits to inner-to-outer order (was outer-to-inner for earlier steps)
    #   in this way, inner clips correspond to junctions, as for merged and aln1
    #   essentially, treat all unmerged split reads the same as a merged split read
    ($alns[0][FLAG] & _SECOND_IN_PAIR) and @alns = reverse(@alns);    
    
    # take junction between each pair of alignments along the read
    foreach my $i2(1..$#alns){ 
        my $i1 = $i2 - 1;
        my @alns = map { [@{$alns[$_]}] } $i1, $i2; # copy to avoid permanently altering underlying alignments by ref             
        setJunction(SPLIT, $alns[READ1], $alns[READ2]);
    }
}
# function for printing junctions in proper order
# here, outer and inner clips are for the _alignments_ immediately flanking the predicted junction
sub setJunction { 
    my ($nodeClass, @alns) = @_;

    # collect endpoint data on the alignments that flank the junction
    # this step assigns node alignment sides based on FF strands
    # from this point forward we use node sides, not FLAG REVERSE, to track orientation
    my (@innData);
    setEndpointData(\@{$innData[READ1]}, $alns[READ1], RIGHT, LEFT);
    setEndpointData(\@{$innData[READ2]}, $alns[READ2], LEFT,  RIGHT);
    
    # get the junction type
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $nodeClass); 
    $jxnType eq PROPER and return;

    # reorient alignment pairs to always reflect left-to-right order on the conjoined genome
    my ($read1, $read2) = (READ1,  READ2);
    if($jxnType eq TRANSLOCATION){
        $alns[$read2][RNAME_INDEX] > $alns[$read1][RNAME_INDEX] or ($read1, $read2) = (READ2, READ1);
    } else {
        abs($innData[$read2][_POS] - $innData[$read1][_POS]) > BIN_SIZE or return;
        $innData[$read2][_POS] > $innData[$read1][_POS] or ($read1, $read2) = (READ2, READ1);
    }

    # print rectified nodes and edges    
    my $node1 = getNodeSignature($alns[$read1], $innData[$read1]);
    my $node2 = getNodeSignature($alns[$read2], $innData[$read2]);
    push @jxns, join("::", $node1, $node2, $jxnType);
}
#===================================================================================================

#===================================================================================================
# utility functions for characterizing SV junctions
#---------------------------------------------------------------------------------------------------
# for these purposes, consider the genome as one contiguous reference unit
#   i.e., with chromosomes conjoined in numerical index order
#       ---chr1---|---chr2---|--- etc.
#   such that ~half of translocations are handled as del/dup, half as inv
#---------------------------------------------------------------------------------------------------
# check to confirm that alignments truly do flank a SV and report its type
sub getJxnType {
    my ($aln1, $aln2, $innData1, $innData2, $nodeClass) = @_;  

    # simplest case, translocation = alignment to different chromosomes
    # PDL vs. inversion sub-handling determined by isCanonicalStrand
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;

    # inversions determined by unexpected alignment orientiations
    $$innData1[_SIDE] eq $$innData2[_SIDE] and return INVERSION;

    # PDL are distinguish based on their position along the conjoined chromosomes
    my $dist = $$innData1[_SIDE] eq LEFTWARD ? 
                    $$innData2[_POS] - $$innData1[_POS] : # LR pair, canonical PDL strand
                    $$innData1[_POS] - $$innData2[_POS];  # RL pair, non-canonical strand (will be flipped later)

    if($nodeClass == GAP){ # based on TLEN=insertSize (since we lack a junction)
        $dist < $gapDupLimit and return DUPLICATION;
        $dist > $gapDelLimit and return DELETION;
    } else { # based on nodes that declare a sequenced junction
        $dist <= 0 and return DUPLICATION;
        $dist > 1  and return DELETION; # pos delta is 1 for a continuous alignment
    }

    # override any aligner designations based on our chrom/strand/separation criteria
    # i.e., we may call some things proper that the aligner did not
    return PROPER;
}
#===================================================================================================

#===================================================================================================
# utility functions for processing individual _alignment_ segments
#---------------------------------------------------------------------------------------------------
# order a set of alignments for a (merged) read in source molecule order
# shortest outer clip is outermost alignment on the read
sub sortReadAlignments {
    my ($alns, $fwdSide, $revSide) = @_;
    my @outDatas;
    foreach my $i(0..$#$alns){
        setEndpointData(\@{$outDatas[$i]}, $$alns[$i], $fwdSide, $revSide);
    }
    return (\@outDatas, sort { $outDatas[$a][_CLIP] <=> $outDatas[$b][_CLIP] } 0..$#$alns);  
}
# collect information on one endpoint (outer or inner) of an alignment
sub setEndpointData {
    my ($data, $aln, $fwdSide, $revSide) = @_;
    my ($clip, $getEnd, $side) = ($$aln[FLAG] & _REVERSE) ?
        ($clips[$revSide], $revSide, $sides[$revSide]) :
        ($clips[$fwdSide], $fwdSide, $sides[$fwdSide]);
    $$data[_CLIP] = $$aln[CIGAR] =~ m/$clip/ ? $1 : 0; # the number of clipped bases
    $$data[_POS]  = $getEnd ? getEnd($$aln[POS], $$aln[CIGAR]) : $$aln[POS]; # the coordinate of the last aligned base
    $$data[_SIDE] = $side; # the direction the aligned read moves away from that coordinate
    $$data[_NODE] = join(":", $$aln[RNAME_INDEX], @$data[_SIDE, _POS]); # complete signature of the SV breakpoint position
}
#===================================================================================================

#===================================================================================================
# collect nodes that identify the outer endpoints of every source molecule
# tally information to create a comprehensive genome fragment coverage map
#---------------------------------------------------------------------------------------------------
sub setOuterEndpoints { 
    my ($outAln1, $outAln2, $outData1, $outData2) = @_; 
    my $node1 = getNodeSignature($outAln1, $outData1);
    my $node2 = getNodeSignature($outAln2, $outData2);
    $molKey = join("::", sort ($node1, $node2));
}
sub getNodeSignature {
    my ($aln, $node) = @_;
    my $pos = $$node[_POS] + ($$node[_SIDE] eq LEFTWARD ? 1 : -1) * $$node[_CLIP]; # hypothetical molecule endpoint at this alignment
    join(":", # this endpoint's signature 
        $$aln[RNAME_INDEX], 
        $$node[_SIDE], 
        $pos < 1 ? 1 : $pos
    );
}
#===================================================================================================

1;
