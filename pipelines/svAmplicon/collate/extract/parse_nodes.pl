use strict;
use warnings;

# constants
use constant {
    AMP_CHROM1 => 0, # amplicon fields
    AMP_SIDE1 => 1,
    AMP_POS1 => 2,
    AMP_CHROM2 => 3, # amplicon fields
    AMP_SIDE2 => 4,
    AMP_POS2 => 5,
    AMP_MOL_COUNT => 6,
    AMP_AMPLICON_ID => 7,
    AMP_INDEX => 8,
    AMP_PROPER => 9,
    #-------------
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    ALN_N => 10,
    RNAME_INDEX => 11,
    #-------------
    _IS_PAIRED => 1, # SAM FLAG bits
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
    #-------------
    MOL_ID => 0,    # molecule-level data, carried here in QNAME
    AMPLICON_ID => 1, 
    N_OVERLAP_BASES => 2, 
    MOL_COUNT => 3, 
    IS_MERGED => 4, 
    READ_N => 5, 
    MOL_CLASS => 5, # values added by extract_nodes regardless of pipeline or bam source
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
    _SHIFT_REVERSE => 4, # how far to shift those bits to yield binary values
    #-------------
    _CLIP => 0, # outData and innData fields
    _POS  => 1,
    _SIDE => 2,
    _SEQ  => 3,
    _NODE => 4,
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
    INSERTION     => 'N', 
    MERGE_FAILURE => 'M',
    #-------------
    _NULL => '*',   
};

# operating parameters
use vars qw($leftClip_ $rightClip_);
my @clips = ($leftClip_, $rightClip_);
my @sides = (RIGHTWARD,  LEFTWARD); # e.g., at a left end clip, the read continues righward from that point

# working variables
use vars qw($READ_LEN $MAX_INSERT_SIZE @alns @mol $jxnN $amplicon);    

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
# source molecules that had one contiguous alignment, only possible in expectOverlap mode
sub commitContiguousAlignment { # aln1/umi1 are the left/proximal end of the source molecule as provided to the aligner
    my ($aln) = @_;
    $mol[MOL_CLASS] = IS_PROPER; # molecule may still have a small indel that did not cause alignment splitting
    setEndpointData(\my @outData1, $aln, LEFT,  RIGHT); # LEFT = process the left clip if mapped on the forward strand
    setEndpointData(\my @outData2, $aln, RIGHT, LEFT);  
    summarizeMolecule([$aln], \@outData1, \@outData2, PROPER); # molecule has no inner endpoints
}
#---------------------------------------------------------------------------------------------------
# a merged read-pair (or single read) aligned in at least two segments, i.e., that crossed an SV junction
sub parseMergedSplit {
    # shortest clip1 will match umi1, longest will match umi2
    my ($outData1s, @alnIs) = sortReadAlignments(\@alns, LEFT, RIGHT);
    $mol[MOL_CLASS] = IS_SV;
    setEndpointData(\my @outData2, $alns[$alnIs[$#alnIs]], RIGHT, LEFT);
    my $jxnTypes = printSequencedJunctions(@alns[@alnIs]); 
    summarizeMolecule([@alns[@alnIs]], $$outData1s[$alnIs[0]], \@outData2, $jxnTypes); # molecule has no inner endpoints
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read pair with no split alignments (i.e., exactly one aln each read)
# either a merge failure or anomalous due to an unaligned junction (possibly at an inner clip)
sub parseUnmergedHiddenJunction {
    my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
    $mol[MOL_CLASS] = IS_SV;
    setEndpointData(\my @innData1, $alns[$read1], RIGHT, LEFT);
    setEndpointData(\my @innData2, $alns[$read2], LEFT,  RIGHT);
    setEndpointData(\my @outData1, $alns[$read1], LEFT,  RIGHT);
    setEndpointData(\my @outData2, $alns[$read2], RIGHT, LEFT);  
    my $jxnType = getJxnType($alns[$read1], $alns[$read2], \@innData1, \@innData2, GAP);
    if($jxnType eq PROPER or $jxnType eq MERGE_FAILURE){
        $mol[MOL_CLASS] = IS_PROPER; 
    } else {   
        # printJunction(GAP, $alns[$read1], $alns[$read2]);     
    }
    summarizeMolecule([$alns[$read1], $alns[$read2]], \@outData1, \@outData2, $jxnType);  
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read-pair that sequenced an SV junction in >=1 of its reads (the most complex alignment patterns land here)
sub parseUnmergedSplit {
    # determine which reads were split
    my (@alnsByRead, @alnIs, @outDatas);
    foreach my $aln(@alns){
        my $read = $$aln[FLAG] & _FIRST_IN_PAIR ? READ1: READ2;
        push @{$alnsByRead[$read]}, $aln;
    }     
    # catch complex pairs where one unmerged read didn't get to us
    ($alnsByRead[READ1] and $alnsByRead[READ2]) or return;
    # catch prior merge failures at over-running read pairs
    if(@{$alnsByRead[READ1]} == @{$alnsByRead[READ2]}){
        my $key1 = join(":", sort map { join(":", @$_[RNAME_INDEX,POS]) } @{$alnsByRead[READ1]});
        my $key2 = join(":", sort map { join(":", @$_[RNAME_INDEX,POS]) } @{$alnsByRead[READ2]});
        if($key1 eq $key2){ # same set of split alignments for each read in an unmerged pair
            @alns = @{$alnsByRead[READ1]}; # eliminate one read, it is identical to the other
            return;
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
    my $outI1 = $alnIs[READ1][0];
    my $outI2 = $alnIs[READ2][0]; 
    $mol[MOL_CLASS] = IS_SV;      
    my $innI1 = $alnIs[READ1][$#{$alnIs[READ1]}];
    my $innI2 = $alnIs[READ2][$#{$alnIs[READ2]}]; 
    my @jxnTypes = printJunction(GAP, $alnsByRead[READ1][$innI1], $alnsByRead[READ2][$innI2]); 
    foreach my $read(READ1, READ2){
        push @jxnTypes, @{$alnsByRead[$read]} > 1 ?
            printSequencedJunctions(@{$alnsByRead[$read]}[@{$alnIs[$read]}]) :
            PROPER;
    }  
    summarizeMolecule(
        [@{$alnsByRead[READ1]}[@{$alnIs[READ1]}], @{$alnsByRead[READ2]}[@{$alnIs[READ2]}]], 
        $outDatas[READ1][$outI1], $outDatas[READ2][$outI2], 
        join("::", @jxnTypes)
    );  
} 
#===================================================================================================

#===================================================================================================
# handle SV junction output
#---------------------------------------------------------------------------------------------------
# a read definitively crosses a junction via one or more split alignments
sub printSequencedJunctions{ 
    my (@alns) = @_; # two or more alignments for a single contiguous (perhaps merged) read
    
    # reorder aln2 splits to inner-to-outer order (was outer-to-inner for earlier steps)
    #   in this way, inner clips correspond to junctions, as for merged and aln1
    #   essentially, treat all unmerged split reads the same as a merged split read
    ($alns[0][FLAG] & _SECOND_IN_PAIR) and @alns = reverse(@alns);    
    
    # take junction between each pair of alignments along the read
    my @jxnTypes;
    foreach my $i2(1..$#alns){ 
        my $i1 = $i2 - 1;

        # copy to avoid permanently altering underlying alignments by ref
        my @alns = map { [@{$alns[$_]}] } $i1, $i2;              
        
        # override read #s to reflect the order of the alignments in the split "pair"
        $alns[READ1][FLAG] |=  _FIRST_IN_PAIR;
        $alns[READ1][FLAG] &= ~_SECOND_IN_PAIR;
        $alns[READ2][FLAG] &= ~_FIRST_IN_PAIR;
        $alns[READ2][FLAG] |=  _SECOND_IN_PAIR;

        # commit the junction
        push @jxnTypes, printJunction(SPLIT, $alns[READ1], $alns[READ2]);
    }
    join(":", @jxnTypes);
}
# function for printing junctions in proper order
# here, outer and inner clips are for the _alignments_ immediately flanking the predicted junction
sub printJunction { 
    my ($nodeClass, @alns) = @_;

    # collect endpoint data on the alignments that flank the junction
    # this step assigns node alignment sides based on FF strands
    # from this point forward we use node sides, not FLAG REVERSE, to track orientation
    my (@innData);
    setEndpointData(\@{$innData[READ1]}, $alns[READ1], RIGHT, LEFT);
    setEndpointData(\@{$innData[READ2]}, $alns[READ2], LEFT,  RIGHT);
    
    # get the junction type
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $nodeClass);  

    return $jxnType;

    # my $isCanonical = isCanonicalStrand($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $jxnType); 
    # $jxnN++;

    # # reorient alignment pairs to always reflect left-to-right order on the canonical strand
    # my ($read1, $read2) = (READ1,  READ2);
    # if(!$isCanonical){
    #     ($read1, $read2) = (READ2,  READ1);
    #     $alns[$read1][FLAG] ^= (_REVERSE + _FIRST_IN_PAIR + _SECOND_IN_PAIR); # for completeness
    #     $alns[$read2][FLAG] ^= (_REVERSE + _FIRST_IN_PAIR + _SECOND_IN_PAIR); # downstream svx tools do not use FLAG
    # }

    # # undo the RC action that the aligner did to the read segment inside an inversion (not the one flanking it)
    # # thus, both SEQs always exit relative to the source molecule on the canonical strand
    # if ($innData[$read1][_SIDE] eq $innData[$read2][_SIDE]){
    #     my $rcRead = ($innData[$read1][_SIDE] eq LEFTWARD) ? $read2 : $read1;
    #     rc(\$alns[$rcRead][SEQ]);
    #     my @cigar = ($alns[$rcRead][CIGAR] =~ m/(\d+\D)/g); # reverse CIGAR also
    #     $alns[$rcRead][CIGAR] = join("", reverse(@cigar));  # FLAG for strand no longer informative!
    # } 

    # # print rectified nodes and edges    
    # printNode($alns[$read1], $innData[$read1], $nodeClass, $jxnType, $jxnN);
    # printNode($alns[$read2], $innData[$read2], $nodeClass, $jxnType, $jxnN);
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
        if($$amplicon[AMP_PROPER] eq "expectOverlap"){ # read pair should have merged, why didn't it?
            $dist <= 0 and return MERGE_FAILURE; # these should be mostly low quality reads
            $dist  > 0 and return INSERTION;     # but novel sequence in the span could also lead failed merging
        } elsif($$amplicon[AMP_PROPER] eq "expectGaps") { # POOR DESIGN CHOICE, ignore until someone needs this path
            return PROPER;
            # my $refSpan = $$amplicon[AMP_POS2] - $$amplicon[AMP_POS1] + 1; # this code is a start but not adequate
            # my $expectedDist = $refSpan - 2 * $READ_LEN + 1;
            # my $gapDelLimit = $MAX_INSERT_SIZE - 2 * $READ_LEN;
            # $dist < $expectedDist and return DELETION;
            # $dist > $expectedDist and return INSERTION;
        } else { # any read pairs from "notPossible" amplicons are SVs of the type consistent with the amplicon design
            return $$amplicon[AMP_POS2] > $$amplicon[AMP_POS1] ? DELETION : DUPLICATION;
        }
    } else { # based on nodes that declare a sequenced junction
        $dist <= 0 and return DUPLICATION;
        $dist  > 1 and return DELETION; # pos delta is 1 for a continuous proper alignment
    }

    # override any aligner designations based on our chrom/strand/separation criteria
    # i.e., we may call some things proper that the aligner did not
    return PROPER;
}
# determine if an alignment pair is from the canonical strand of a junction fragment
# this is a junction-level property and thus could differ from the moleculeStrand
#   the canonical strand is:
#       the top genome strand for PDL and tranlocations handled as PDL
#       the top genome strand for the alignment NOT in the inverted segment
#           e.g., the top strand to the left of the left inversion junction and vice versa
sub isCanonicalStrand {
    my ($aln1, $aln2, $innData1, $innData2, $jxnType) = @_;
    if($$innData1[_SIDE] eq $$innData2[_SIDE]){ # inversion-type handling, sort by pos along conjoined chromosomes
        if($jxnType eq INVERSION){ # inversions
            $$innData1[_POS] < $$innData2[_POS];
        } else { # inversion-type translocations
            $$aln1[RNAME_INDEX] < $$aln2[RNAME_INDEX];
        }
    } else { # del/dup and PDL-type translocations, canonical is LR pairs
        $$innData1[_SIDE] eq LEFTWARD;
    }
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
    $$data[_SEQ] = $$data[_CLIP] > 0 ? # the sequence of the clipped bases
        ($getEnd ?
            substr($$aln[SEQ], -$$data[_CLIP]) :
            substr($$aln[SEQ], 0, $$data[_CLIP])):
        _NULL;
    $$data[_NODE] = join(":", $$aln[RNAME_INDEX], @$data[_SIDE, _POS]); # complete signature of the SV breakpoint position
}
#===================================================================================================

#===================================================================================================
# print molecule-defining data bits
#---------------------------------------------------------------------------------------------------
sub summarizeMolecule {
    my ($orderedAlns, $outData1, $outData2, $jxnType) = @_; 
    my $nBases = $mol[IS_MERGED] ? 
        length($$orderedAlns[0][SEQ]) :
        $mol[N_OVERLAP_BASES] eq "NA" ?
            "NA" :
            $READ_LEN + $READ_LEN - $mol[N_OVERLAP_BASES];
    print join("\t",
        @mol, 
        $jxnType,
        @$outData1[_POS], #
        @$outData2[_POS], 
        scalar(@$orderedAlns),
        $nBases,
        join("::", map { join(":", @$_[FLAG, RNAME, POS,MAPQ, CIGAR]) } @$orderedAlns)
    ), "\n"; 
}
#===================================================================================================

#===================================================================================================
# print candidate SV data per source molecule:
#   nodes = alignment endpoints at alignment discontinuities (clips, splits or gaps)
#   edges = SV junctions that connect two nodes (whether sequenced or just inferred)
#---------------------------------------------------------------------------------------------------
# all available SV evidence nodes, listed one row per node over all source molecules
sub printNode { 
    my ($aln, $node, $nodeClass, $jxnType, $jxnN) = @_;
    my $alnN = $mol[MOL_CLASS] eq IS_PROPER ? 0 : $$aln[ALN_N]; # even unmerged proper molecules connect their outer nodes
    # $nodesH
    # print join("\t", # nodeN (1,2) can be inferred from order in the nodes file
    #     #------------------------------------------- # below this line is (potentially) unique to each node
    #     @$node[_NODE, _CLIP, _SEQ],                  # node-level data
    #     @$aln[FLAG, POS, MAPQ, CIGAR, SEQ], $alnN,   # alignment-level data
    #     #------------------------------------------- # below this line is common to both nodes in a junction
    #     $nodeClass,                                  # node-level data
    #     $jxnType, $jxnN,                             # edge/junction-level data 
    #     @mol[MOL_ID..MOL_CLASS]                      # molecule-level data
    # ), "\n";
}
#===================================================================================================

1;
