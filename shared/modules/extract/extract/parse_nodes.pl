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
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    ALN_N => 10,
    RNAME_INDEX => 11,
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
    MAX_FAMILY_COUNT => 100
};

# operating parameters
use vars qw($leftClip_ $rightClip_);
my @clips = ($leftClip_, $rightClip_);
my @sides = (RIGHTWARD,  LEFTWARD); # e.g., at a left end clip, the read continues righward from that point

# working variables
use vars qw(@alns @mol $jxnN
            $spansH $nodesH $endpointsH
            $isTargeted $maxInsertSize %insertSizes $isCountStrands %strandCounts 
            $fwdSide2 $revSide2
            $IS_COLLATED $LIBRARY_TYPE $MIN_CLIP $MAX_TLEN $READ_LEN);    
my $gapDelLimit =  2 * $MAX_TLEN - 2 * $READ_LEN; # these are conservative, i.e., call more proper than aligner typically
my $gapDupLimit = -2 * $READ_LEN;
my $collapseStrands = ($LIBRARY_TYPE eq 'TruSeq'); # otherwise, Nextera, where same signatures on opposite strands are unique source molecules
my $isFR = !$IS_COLLATED;
my (@outerPos, @outerPosWrk);

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
# source molecules that are considered proper (merged or unmerged)
sub commitProperMolecule { # aln1/umi1 are the left/proximal end of the source molecule as provided to the aligner
    my ($aln1, $aln2, $fwdSide2, $revSide2) = @_;
    $mol[MOL_CLASS] = IS_PROPER;
    !defined $mol[MOL_STRAND] and $mol[MOL_STRAND] = ($$aln1[FLAG] & _REVERSE) >> _SHIFT_REVERSE; # 0=forward, 1=reverse
    setEndpointData(\my @outData1, $aln1, LEFT,  RIGHT); # LEFT = process the left clip if mapped on the forward strand
    setEndpointData(\my @outData2, $aln2, $fwdSide2, $revSide2); # orientation depends on merge state
    $mol[IS_OUTER_CLIP1] = ($outData1[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outData2[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $$aln1[RNAME],   $$aln2[RNAME], 
        $outData1[_POS], $outData2[_POS]
    );
    printOuterEndpoints($aln1, $aln2, \@outData1, \@outData2);  
    # prepare to write a cross-tabulated file of proper molecule insert sizes
    my $insertSize = abs($outData2[_POS] - $outData1[_POS]) + 1 + $outData1[_CLIP] + $outData2[_CLIP];
    $insertSize <= $maxInsertSize and $insertSizes{$mol[TARGET_CLASS]}[$insertSize]++;
}
#---------------------------------------------------------------------------------------------------
# a merged read-pair (or single read) aligned in at least two segments, i.e., that crossed an SV junction
sub parseMergedSplit {
    # shortest clip1 will match umi1, longest will match umi2
    my ($outData1s, @alnIs) = sortReadAlignments(\@alns, LEFT, RIGHT);
    $mol[MOL_CLASS] = IS_SV;
    !defined $mol[MOL_STRAND] and $mol[MOL_STRAND] = getMoleculeStrand(@alns[$alnIs[0], $alnIs[$#alnIs]]);
    setEndpointData(\my @outData2, $alns[$alnIs[$#alnIs]], RIGHT, LEFT);
    $mol[IS_OUTER_CLIP1] = ($$outData1s[$alnIs[0]][_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outData2[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $alns[$alnIs[0]][RNAME],      $alns[$alnIs[$#alnIs]][RNAME], 
        $$outData1s[$alnIs[0]][_POS], $outData2[_POS]
    );
    printOuterEndpoints(@alns[$alnIs[0], $alnIs[$#alnIs]], $$outData1s[$alnIs[0]], \@outData2);     
    printSequencedJunctions(@alns[@alnIs]); 
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read pair with no split alignments (i.e., exactly one aln each read)
# but anomalous due to an unaligned junction (possibly at an inner clip)
sub parseUnmergedHiddenJunction {
    my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
    $mol[MOL_CLASS] = IS_SV;
    !defined $mol[MOL_STRAND] and $mol[MOL_STRAND] = getMoleculeStrand($alns[$read1], $alns[$read2]);
    # override aligner proper flag sometimes (aligner fails to call some things proper, e.g., failed-merge overruns)
    setEndpointData(\my @innData1, $alns[$read1], RIGHT, LEFT);
    setEndpointData(\my @innData2, $alns[$read2], $revSide2, $fwdSide2);
    my $jxnType = getJxnType($alns[$read1], $alns[$read2], \@innData1, \@innData2, GAP);
    $jxnType eq PROPER and return commitProperMolecule($alns[$read1], $alns[$read2], $fwdSide2, $revSide2);    
    # commit SV if truly not proper    
    setEndpointData(\my @outData1, $alns[$read1], LEFT,  RIGHT);
    setEndpointData(\my @outData2, $alns[$read2], $fwdSide2, $revSide2);
    $mol[IS_OUTER_CLIP1] = ($outData1[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outData2[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $alns[$read1][RNAME], $alns[$read2][RNAME], 
        $outData1[_POS],      $outData2[_POS]
    );
    printOuterEndpoints($alns[$read1], $alns[$read2], \@outData1, \@outData2); # actions based on inner portion of molecule    
    printJunction(GAP, $revSide2, $fwdSide2, $alns[$read1], $alns[$read2]);      
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
        my ($fwdSide, $revSide) = ($read == READ1) ? (LEFT, RIGHT) : ($fwdSide2, $revSide2);
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
    $mol[MOL_CLASS] = IS_SV;     
    !defined $mol[MOL_STRAND] and $mol[MOL_STRAND] = getMoleculeStrand($alnsByRead[READ1][$outI1], $alnsByRead[READ2][$outI2]);
    $mol[IS_OUTER_CLIP1] = ($outDatas[READ1][$outI1][_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outDatas[READ2][$outI2][_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $alnsByRead[READ1][$outI1][RNAME], $alnsByRead[READ2][$outI2][RNAME], 
        $outDatas[READ1][$outI1][_POS],    $outDatas[READ2][$outI2][_POS]
    );
    printOuterEndpoints(
        $alnsByRead[READ1][$outI1], $alnsByRead[READ2][$outI2],
          $outDatas[READ1][$outI1],   $outDatas[READ2][$outI2]
    );    
    # actions based on inner portion of molecule
    # all junctions in all split reads
    foreach my $read(READ1, READ2){
        @{$alnsByRead[$read]} > 1 and
            printSequencedJunctions(@{$alnsByRead[$read]}[@{$alnIs[$read]}]);
    }  
    # possible additional SV junction in the gap
    my $innI1 = $alnIs[READ1][$#{$alnIs[READ1]}];
    my $innI2 = $alnIs[READ2][$#{$alnIs[READ2]}]; 
    printJunction(GAP, $revSide2, $fwdSide2, $alnsByRead[READ1][$innI1], $alnsByRead[READ2][$innI2]); 
} 
#===================================================================================================

#===================================================================================================
# handle SV junction output
#---------------------------------------------------------------------------------------------------
# a read definitively crosses a junction via one or more split reads
sub printSequencedJunctions{ 
    my (@alns) = @_; # two or more alignments for a single contiguous (perhaps merged) read
    
    # reorder aln2 splits to inner-to-outer order (was outer-to-inner for earlier steps)
    #   in this way, inner clips correspond to junctions, as for merged and aln1
    #   essentially, treat all unmerged split reads the same as a merged split read
    ($alns[0][FLAG] & _SECOND_IN_PAIR) and @alns = reverse(@alns);    
    
    # take junction between each pair of alignments along the read
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
        printJunction(SPLIT, LEFT, RIGHT, $alns[READ1], $alns[READ2]);
    }
}
# function for printing junctions in proper order
# here, outer and inner clips are for the _alignments_ immediately flanking the predicted junction
sub printJunction { 
    my ($nodeClass, $innFwdSide2, $innRevSide2, @alns) = @_;

    # collect endpoint data on the alignments that flank the junction
    my (@innData);
    setEndpointData(\@{$innData[READ1]}, $alns[READ1], RIGHT, LEFT);
    setEndpointData(\@{$innData[READ2]}, $alns[READ2], $innFwdSide2, $innRevSide2);
    
    # get the junction type
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $nodeClass);  
    my $isCanonical = isCanonicalStrand($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $jxnType); 
    $jxnN++;

    # reorient alignment pairs to always reflect left-to-right order on the canonical strand
    my ($read1, $read2, $umi1, $umi2, $isOuterClip1,   $isOuterClip2) =
        (READ1,  READ2,  UMI1,  UMI2,  IS_OUTER_CLIP1, IS_OUTER_CLIP2);
    if($isCanonical){
        @outerPosWrk = @outerPos;
    } else {
        @outerPosWrk = reverse @outerPos;
        ($read1, $read2, $umi1, $umi2, $isOuterClip1,  $isOuterClip2) =
         (READ2,  READ1,  UMI2,  UMI1,  IS_OUTER_CLIP2, IS_OUTER_CLIP1);
        $alns[$read1][FLAG] ^= (_FIRST_IN_PAIR + _SECOND_IN_PAIR);
        $alns[$read2][FLAG] ^= (_FIRST_IN_PAIR + _SECOND_IN_PAIR); 
        unless($isFR and $nodeClass == GAP){ # downstream SVX tools don't rely on FLAG
            $alns[$read1][FLAG] ^= (_REVERSE);
            $alns[$read2][FLAG] ^= (_REVERSE); 
        }
    }

    # undo the RC action that the aligner did to the read segment inside an inversion (not the one flanking it)
    # thus, both SEQs always exit relative to the source molecule on the canonical strand
    if ($innData[$read1][_SIDE] eq $innData[$read2][_SIDE]){
        my $rcRead = ($innData[$read1][_SIDE] eq LEFTWARD) ? $read2 : $read1;
        rc(\$alns[$rcRead][SEQ]);
        my @cigar = ($alns[$rcRead][CIGAR] =~ m/(\d+\D)/g); # reverse CIGAR also
        $alns[$rcRead][CIGAR] = join("", reverse(@cigar));  # FLAG for strand no longer informative!
    } 

    # print rectified nodes and edges    
    printNode($umi1, $alns[$read1], $innData[$read1], 
              $nodeClass, $jxnType, $jxnN);
    printNode($umi2, $alns[$read2], $innData[$read2], 
              $nodeClass, $jxnType, $jxnN);
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
# assign a source strand to help track uncollated TruSeq strand duplicates of the same source molecule
# this is a molecule-level property calculated on the outermost, i.e., reference alignments
# if previously collated, MOL_STRAND was already set upstream to 0, 1, or 2 (duplex, i.e., both strands)
sub getMoleculeStrand {
    my ($outAln1, $outAln2) = @_;
    my $strand1 = ($$outAln1[FLAG] & _REVERSE) >> _SHIFT_REVERSE; # 0=forward, 1=reverse
    my $strand2 = ($$outAln2[FLAG] & _REVERSE) >> _SHIFT_REVERSE;
    $$outAln2[FLAG] & _IS_PAIRED and $strand2 = ($strand2 + 1) % 2;
    if($strand1 == $strand2){ # PDL and some T
        return $strand1; 
    } elsif($$outAln1[RNAME_INDEX] != $$outAln2[RNAME_INDEX]){ # T, sort by chrom
        return $$outAln1[RNAME_INDEX] > $$outAln2[RNAME_INDEX] ? 1 : 0; 
    } else { # I, sort by pos
        return $$outAln1[POS] > $$outAln2[POS] ? 1 : 0; 
    }  
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
    $$data[_SEQ] = $$data[_CLIP] >= $MIN_CLIP ? # the sequence of the clipped bases
        ($getEnd ?
            substr($$aln[SEQ], -$$data[_CLIP]) :
            substr($$aln[SEQ], 0, $$data[_CLIP])):
        _NULL;
    $$data[_NODE] = join(":", $$aln[RNAME_INDEX], @$data[_SIDE, _POS]); # complete signature of the SV breakpoint position
}
#===================================================================================================

#===================================================================================================
# collect nodes that identify the outer endpoints of every source molecule
# tally information to create a comprehensive genome fragment coverage map
#---------------------------------------------------------------------------------------------------
sub printOuterEndpoints { 
    my ($outAln1, $outAln2, $outData1, $outData2) = @_; 

    # generate two identifying outer positions per source molecule
    # these are what the node position would have been had any outer clip been aligned
    @outerPosWrk = @outerPos = (
        $$outData1[_POS] + ($$outData1[_SIDE] eq 'L' ? 1 : -1) * $$outData1[_CLIP],
        $$outData2[_POS] + ($$outData2[_SIDE] eq 'L' ? 1 : -1) * $$outData2[_CLIP],
    );

    # prepare for setting SHARED_PROPER downstream
    my $node1 = getNodeSignature($outAln1, UMI1, $outData1);
    my $node2 = getNodeSignature($outAln2, UMI2, $outData2);    
    if($endpointsH){
        my $mol = join("\t", @mol[MOL_ID, MOL_CLASS]); # molecule properties
        print $endpointsH join("\t", $node1, $mol), "\n";
        print $endpointsH join("\t", $node2, $mol), "\n";
    }

    # prepare for SV evidence support via clipped nodes
    my $molIsOuterClipped = ($mol[IS_OUTER_CLIP1] or $mol[IS_OUTER_CLIP2]);
    printNode(UMI1, $outAln1, $outData1, # outer clip output used as SV evidence support
              OUTER_CLIP, _NULL, -1, $molIsOuterClipped);
    printNode(UMI2, $outAln2, $outData2, # outer clips bear negative JXN_N
              OUTER_CLIP, _NULL, -1, $molIsOuterClipped);

    # prepare for coverage map analysis downstream if untargeted
    $spansH and print $spansH join("\t", 
        $mol[MOL_CLASS],
        $collapseStrands ? 0 : $mol[MOL_STRAND], # molecule strand (not read-pair strand) and outer node signatures used for read-pair de-duplication
        $mol[MOL_STRAND] ? $node2 : $node1, # flip order for bottom strand proper
        $mol[MOL_STRAND] ? $node1 : $node2,
        $mol[MOL_CLASS] eq IS_PROPER ? 
            ($molIsOuterClipped ? 
                join(":", $mol[MOL_STRAND] ? ($$outData2[_POS], $$outData1[_POS]) : ($$outData1[_POS], $$outData2[_POS])) : 
                "-") : 
            join("::", map {join(":", $$_[RNAME_INDEX], $$_[POS] - 1, getEnd($$_[POS], $$_[CIGAR])) } @alns)
    ), "\n";

    # prepare data to write a cross-tabulated file of _molecule_ strand1 count by strand2 count if targeted
    # stratified by where the _outer_ molecule ends aligned, and if they were proper or SV
    $isCountStrands and $strandCounts{"$mol[TARGET_CLASS]\t$mol[MOL_CLASS]"}[ # strand count order does not matter for crosstab
        $mol[STRAND_COUNT1] > MAX_FAMILY_COUNT ? MAX_FAMILY_COUNT : $mol[STRAND_COUNT1]
    ][
        $mol[STRAND_COUNT2] > MAX_FAMILY_COUNT ? MAX_FAMILY_COUNT : $mol[STRAND_COUNT2]
    ]++;
}
sub getNodeSignature {
    my ($aln, $umi, $node) = @_;
    my $pos = $$node[_POS] + ($$node[_SIDE] eq LEFTWARD ? 1 : -1) * $$node[_CLIP]; # hypothetical molecule endpoint at this alignment
    join(":", # this endpoint's signature 
        $mol[$umi], 
        $$aln[RNAME_INDEX], 
        $$node[_SIDE], 
        $pos < 1 ? 1 : $pos, 
        $$node[_CLIP] >= $MIN_CLIP ? 1 : 0  # flag whether a refAln pos was clipped
    );
}
#===================================================================================================

#===================================================================================================
# print candidate SV data per source molecule:
#   nodes = alignment endpoints at alignment discontinuities (clips, splits or gaps)
#   edges = SV junctions that connect two nodes (whether sequenced or just inferred)
#---------------------------------------------------------------------------------------------------
# all available SV evidence nodes, listed one row per node over all source molecules
sub printNode { 
    my ($umi, $aln, $node, 
        $nodeClass, $jxnType, $jxnN, $molIsOuterClipped) = @_;
    if($nodeClass == OUTER_CLIP and # all internal junction nodes are printed
       !$molIsOuterClipped and      # proper molecules with one outer clip are printed at both ends
       $mol[MOL_CLASS] eq IS_PROPER # all nodes on SV molecules are printed, even unclipped outer
    ){ return } # don't print anything for unclipped proper molecules
    my $alnN = $mol[MOL_CLASS] eq IS_PROPER ? 0 : $$aln[ALN_N]; # even unmerged proper molecules connect their outer nodes
    print $nodesH join("\t", # nodeN (1,2) can be inferred from order in the nodes file
        #------------------------------------------- # below this line is (potentially) unique to each node
        @$node[_NODE, _CLIP, _SEQ],                  # node-level data
        @$aln[FLAG, POS, MAPQ, CIGAR, SEQ], $alnN,   # alignment-level data
        $mol[$umi],                                  # molecule-level data
        #------------------------------------------- # below this line is common to both nodes in a junction
        $nodeClass,                                  # node-level data
        $jxnType, $jxnN,                             # edge/junction-level data 
        @mol[MOL_ID, IS_MERGED..SHARED_PROPER],      # molecule-level data
        @outerPosWrk # first OUT_POS matches the first clip node in the file as written
    ), "\n";
}
#===================================================================================================

1;
