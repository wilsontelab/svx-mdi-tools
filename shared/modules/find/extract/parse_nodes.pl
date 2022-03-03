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
            $IS_COLLATED $TARGETS_BED $LIBRARY_TYPE $MIN_CLIP $MAX_TLEN $READ_LEN);    
my $gapDelLimit =  2 * $MAX_TLEN - 2 * $READ_LEN; # these are conservative, i.e., call more proper than aligner typically
my $gapDupLimit = -2 * $READ_LEN;
my $collapseStrands = ($LIBRARY_TYPE eq 'TruSeq'); # otherwise, Nextera, where same signatures on opposite strands are unique source molecules

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
# source molecules that are considered proper (merged or unmerged)
sub commitProperMolecule { # aln1/umi1 are the left/proximal end of the source molecule as provided to the aligner
    my ($aln1, $aln2, $fwdSide2, $revSide2) = @_;
    $mol[MOL_CLASS] = IS_PROPER;
    $mol[MOL_STRAND] = ($$aln1[FLAG] & _REVERSE) >> _SHIFT_REVERSE; # 0=forward, 1=reverse
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
    $mol[MOL_STRAND] = getMoleculeStrand(@alns[$alnIs[0], $alnIs[$#alnIs]]);
    setEndpointData(\my @outData2, $alns[$alnIs[$#alnIs]], RIGHT, LEFT);
    $mol[IS_OUTER_CLIP1] = ($$outData1s[$alnIs[0]][_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outData2[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $alns[$alnIs[0]][RNAME],      $alns[$alnIs[$#alnIs]][RNAME], 
        $$outData1s[$alnIs[0]][_POS], $outData2[_POS]
    );
    printSequencedJunctions(@alns[@alnIs]); 
    printOuterEndpoints(@alns[$alnIs[0], $alnIs[$#alnIs]], $$outData1s[$alnIs[0]], \@outData2);    
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read pair with no split alignments (i.e., exactly one aln each read)
# but anomalous due to an unaligned junction (possibly at an inner clip)
sub parseUnmergedHiddenJunction {
    my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
    $mol[MOL_CLASS] = IS_SV;
    $mol[MOL_STRAND] = getMoleculeStrand($alns[$read1], $alns[$read2]);
    # override aligner proper flag sometimes (aligner fails to call some things proper, e.g., failed-merge overruns)
    setEndpointData(\my @innData1, $alns[$read1], RIGHT, LEFT);
    setEndpointData(\my @innData2, $alns[$read2], RIGHT, LEFT);
    my $jxnType = getJxnType($alns[$read1], $alns[$read2], \@innData1, \@innData2, GAP);
    $jxnType eq PROPER and return commitProperMolecule($alns[$read1], $alns[$read2], LEFT,  RIGHT);    
    # commit SV if truly not proper    
    setEndpointData(\my @outData1, $alns[$read1], LEFT,  RIGHT);
    setEndpointData(\my @outData2, $alns[$read2], LEFT,  RIGHT);
    $mol[IS_OUTER_CLIP1] = ($outData1[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outData2[_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $alns[$read1][RNAME], $alns[$read2][RNAME], 
        $outData1[_POS],      $outData2[_POS]
    );
    printJunction(GAP, RIGHT, LEFT, $alns[$read1], $alns[$read2]);      
    printOuterEndpoints($alns[$read1], $alns[$read2], \@outData1, \@outData2); # actions based on inner portion of molecule
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
        # my ($fwdSide, $revSide) = ($read == READ1) ? (LEFT, RIGHT) : (RIGHT, LEFT);
        if(@{$alnsByRead[$read]} == 1){ # an unsplit read paired with a split read
            @{$alnIs[$read]} = (0); # the index within @alnsByRead, not @alns
            setEndpointData(\my @outData, $alnsByRead[$read][0], LEFT, RIGHT);
            $outDatas[$read] = [\@outData];
        } else { # one of either one or two split reads
            ($outDatas[$read], @{$alnIs[$read]}) = sortReadAlignments($alnsByRead[$read], LEFT, RIGHT); 
        }  
    }
    # actions based on outer portion of total source molecule
    my $outI1 = $alnIs[READ1][0];
    my $outI2 = $alnIs[READ2][0]; 
    $mol[MOL_CLASS] = IS_SV;     
    $mol[MOL_STRAND] = getMoleculeStrand($alnsByRead[READ1][$outI1], $alnsByRead[READ2][$outI2]);
    $mol[IS_OUTER_CLIP1] = ($outDatas[READ1][$outI1][_CLIP] >= $MIN_CLIP ? 1 : 0);
    $mol[IS_OUTER_CLIP2] = ($outDatas[READ2][$outI2][_CLIP] >= $MIN_CLIP ? 1 : 0);
    $isTargeted and $mol[TARGET_CLASS] = getTargetClass(
        $alnsByRead[READ1][$outI1][RNAME], $alnsByRead[READ2][$outI2][RNAME], 
        $outDatas[READ1][$outI1][_POS],    $outDatas[READ2][$outI2][_POS]
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
    printJunction(GAP, RIGHT, LEFT, $alnsByRead[READ1][$innI1], $alnsByRead[READ2][$innI2]); 
    printOuterEndpoints(
        $alnsByRead[READ1][$outI1], $alnsByRead[READ2][$outI2],
          $outDatas[READ1][$outI1],   $outDatas[READ2][$outI2]
    );
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

    # collect endpoint data on the alignments that flank the splitting junction
    my (@innData);
    setEndpointData(\@{$innData[READ1]}, $alns[READ1], RIGHT, LEFT);
    setEndpointData(\@{$innData[READ2]}, $alns[READ2], $innFwdSide2, $innRevSide2);
    
    # get the junction type
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $nodeClass);  
    $jxnN++;

    # reorient alignment pairs to always reflect the canonical junction strand
    my ($read1, $read2, $isOuterClip1, $isOuterClip2) =
        (READ1, READ2, IS_OUTER_CLIP1, IS_OUTER_CLIP2);
    my $isCanonical = isCanonicalStrand($jxnType, $nodeClass, @alns);   
    if(!$isCanonical){
        ($read1, $read2, $isOuterClip1, $isOuterClip2) =
        (READ2, READ1, IS_OUTER_CLIP2, IS_OUTER_CLIP1);
        $alns[$read1][FLAG] ^= (_REVERSE + _FIRST_IN_PAIR + _SECOND_IN_PAIR); # TODO: handle _REVERSE differently?? *****
        $alns[$read2][FLAG] ^= (_REVERSE + _FIRST_IN_PAIR + _SECOND_IN_PAIR);    
    }    
    
    # undo the RC action that aligner did to one of the pair of our inversion alignments
    # thus, both SEQs always exit relative to the source molecule on the canonical strand
    if($nodeClass == GAP ?
       ($alns[$read1][FLAG] & _REVERSE) == ($alns[$read2][FLAG] & _REVERSE) :
       ($alns[$read1][FLAG] & _REVERSE) != ($alns[$read2][FLAG] & _REVERSE)){
        my $rcRead = $isCanonical ? # TODO: this may still need fixing for gaps?
            (($alns[$read1][FLAG] & _REVERSE) ? $read1 : $read2) :
            (($alns[$read1][FLAG] & _REVERSE) ? $read2 : $read1);
        rc(\$alns[$rcRead][SEQ]);
        my @cigar = ($alns[$rcRead][CIGAR] =~ m/(\d+\D)/g); # reverse CIGAR also
        $alns[$rcRead][CIGAR] = join("", reverse(@cigar));  # FLAG for strand no longer informative!
    } 

    # print rectified nodes and edges    
    printNode($alns[$read1], $innData[$read1], $nodeClass, $jxnType, $jxnN);
    printNode($alns[$read2], $innData[$read2], $nodeClass, $jxnType, $jxnN);
}
#===================================================================================================

#===================================================================================================
# utility functions for characterizing SV junctions
#---------------------------------------------------------------------------------------------------
# check to confirm that alignments truly do flank a SV and report its type
sub getJxnType {
    my ($aln1, $aln2, $innData1, $innData2, $nodeClass) = @_;  
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;
    my $rev1 = $$aln1[FLAG] & _REVERSE;
    if($nodeClass == GAP){ # based on TLEN=insertSize (since we lack a junction); read strands inverted relative to each other in FR pair
        $rev1 == ($$aln2[FLAG] & _REVERSE) and return INVERSION; 
        my $dist = $rev1 ? ($$innData1[_POS] - $$innData2[_POS]) : ($$innData2[_POS] - $$innData1[_POS]);
        $dist < $gapDupLimit and return DUPLICATION;
        $dist > $gapDelLimit and return DELETION;
    } else { # based on nodes that declare a sequenced junction; reads on the same strand
        $rev1 != ($$aln2[FLAG] & _REVERSE) and return INVERSION; # given that we expect FF or RR after for a sequenced junction
        my $dist = $rev1 ? ($$innData1[_POS] - $$innData2[_POS]) : ($$innData2[_POS] - $$innData1[_POS]);
        $dist <= 0 and return DUPLICATION;
        $dist > 1 and return DELETION;
    }
    return PROPER;
}
# assign a source strand to help track TruSeq strand duplicates of the same source molecule
# this is a molecule-level property calculated on the outermost, i.e., reference alignments
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
sub isCanonicalStrand {
    my ($jxnType, $nodeClass, $aln1, $aln2) = @_;  
    if($jxnType eq DELETION or 
       $jxnType eq DUPLICATION){
        return $$aln1[FLAG] & _REVERSE ? 0 : 1;
    } elsif($jxnType eq INVERSION) {
        return $$aln1[POS] > $$aln2[POS] ? 0 : 1;
    } elsif($jxnType eq TRANSLOCATION){
        if($nodeClass == GAP ? 
           ($$aln1[FLAG] & _REVERSE) != ($$aln2[FLAG] & _REVERSE) : 
           ($$aln1[FLAG] & _REVERSE) == ($$aln2[FLAG] & _REVERSE)){ 
            return $$aln1[FLAG] & _REVERSE ? 0 : 1;
        } else {
            return $$aln1[RNAME_INDEX] > $$aln2[RNAME_INDEX] ? 0 : 1;
        }
    } else { # unknown, extremely rare
        return 1;
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
    $$data[_NODE] = join(":", $$aln[RNAME_INDEX], @$data[_SIDE, _POS]); # signature of the SV breakpoint position
}
#===================================================================================================

#===================================================================================================
# collect nodes that identify the outer endpoints of every source molecule
# tally information to create a comprehensive genome fragment coverage map
#---------------------------------------------------------------------------------------------------
sub printOuterEndpoints { 
    my ($outAln1, $outAln2, $outData1, $outData2) = @_; 

    # prepare for setting SHARED_PROPER downstream
    if($endpointsH){
        printEndpoint($outAln1, UMI1, $outData1); # output used to do proper molecule matching
        printEndpoint($outAln2, UMI2, $outData2);
    }

    # prepare for SV evidence support via clipped nodes
    my $molIsOuterClipped = ($mol[IS_OUTER_CLIP1] or $mol[IS_OUTER_CLIP2]);
    printNode(UMI1, $outAln1, $outData1, OUTER_CLIP, _NULL, -1, $molIsOuterClipped); # outer clip output used as SV evidence support
    printNode(UMI2, $outAln2, $outData2, OUTER_CLIP, _NULL, -2, $molIsOuterClipped); # outer clips bear negative JXN_N, 1/2 for each molecule end

    # prepare for coverage map analysis downstream if untargeted
    $spansH and print $spansH join("\t", 
        $mol[MOL_CLASS],
        $collapseStrands ? 0 : $mol[MOL_STRAND], # molecule strand (not read-pair strand) and outer node signatures used for read-pair de-duplication
        $mol[MOL_STRAND] ? getNodeSignature($outAln2, $outData2) : getNodeSignature($outAln1, $outData1), # flip order for bottom strand proper
        $mol[MOL_STRAND] ? getNodeSignature($outAln1, $outData1) : getNodeSignature($outAln2, $outData2),
        $mol[MOL_CLASS] eq IS_PROPER ? 
            ($molIsOuterClipped ? getProperSpan($outData1, $outData2) : "-") : 
            getAlignmentSpans() # information for creating the coverage map
    ), "\n";

    # prepare data to write a cross-tabulated file of _molecule_ strand1 count by strand2 count if targeted
    # stratified by where the _outer_ molecule ends aligned, and if they were proper or SV
    $isCountStrands and $strandCounts{"$mol[TARGET_CLASS]\t$mol[MOL_CLASS]"}[ # strand count order does not matter for crosstab
        $mol[STRAND_COUNT1] > MAX_FAMILY_COUNT ? MAX_FAMILY_COUNT : $mol[STRAND_COUNT1]
    ][
        $mol[STRAND_COUNT2] > MAX_FAMILY_COUNT ? MAX_FAMILY_COUNT : $mol[STRAND_COUNT2]
    ]++;
}
sub printEndpoint { # all source molecule outer endpoint nodes, whether SV evidence or not
    my ($aln, $umi, $node) = @_;
    my $sign = $$node[_SIDE] eq LEFTWARD ? 1 : -1;    
    my $pos = $$node[_POS] + $sign * $$node[_CLIP];        # hypothetical molecule endpoint at this alignment
    my $isSVClip = $$node[_CLIP] >= $MIN_CLIP ? 1 : 0;     # flag whether a refAln pos was clipped
    print $endpointsH join("\t",                           # overall, same endpoint ID strategy as read pair grouping
        @mol[MOL_ID, MOL_CLASS, MOL_STRAND, TARGET_CLASS], # molecule properties
        join(":", $mol[$umi], $$aln[RNAME_INDEX], $$node[_SIDE], $pos, $isSVClip) # this endpoint's signature
    ), "\n";  
}
sub getNodeSignature {
    my ($aln, $node) = @_;
    my $sign = $$node[_SIDE] eq LEFTWARD ? 1 : -1;    
    my $pos = max(1, $$node[_POS] + $sign * $$node[_CLIP]); # hypothetical molecule endpoint at this alignment
    join("\t", $$aln[RNAME_INDEX], $$node[_SIDE], $pos);  # this endpoint's signature
}
sub getProperSpan {
    my ($outData1, $outData2) = @_;
    join(":", 
        $mol[MOL_STRAND] ? 
            ($$outData2[_POS], $$outData1[_POS]) : 
            ($$outData1[_POS], $$outData2[_POS])
    )    
}
sub getAlignmentSpans { # individual BED3 formatted spans for all SV alignments
    join("::", map {
        join(":", $$_[RNAME_INDEX], $$_[POS] - 1, getEnd($$_[POS], $$_[CIGAR])) 
    } @alns);
}
#===================================================================================================

#===================================================================================================
# print candidate SV data per source molecule:
#   nodes = alignment endpoints at alignment discontinuities (clips, splits or gaps)
#   edges = SV junctions that connect two nodes (whether sequenced or just inferred)
#---------------------------------------------------------------------------------------------------
# all available SV evidence nodes, listed one row per node over all source molecules
# carries information for subsequent collapse into junction and alignment edges
sub printNode { 
    my ($umi, $aln, $node, $nodeClass, $jxnType, $jxnN, $molIsOuterClipped) = @_;
    if($nodeClass == OUTER_CLIP and # all internal junction nodes are printed
       !$molIsOuterClipped and      # proper molecules with one outer clip are printed at both ends
       $mol[MOL_CLASS] eq IS_PROPER # all nodes on SV molecules are printed, even unclipped outer
    ){ return } # don't print anything for unclipped proper molecules
    my $alnN = $mol[MOL_CLASS] eq IS_PROPER ? 0 : $$aln[ALN_N]; # even unmerged proper molecules connect their outer nodes
    print $nodesH join("\t", 
        @$node[_NODE, _CLIP, _SEQ], $nodeClass,      # node-level data
        $jxnType, $jxnN,                             # edge/junction-level data
        @$aln[FLAG, POS, MAPQ, CIGAR, SEQ], $alnN,   # alignment-level data
        @mol[MOL_ID, $umi, IS_MERGED..SHARED_PROPER] # molecule-level data
    ), "\n";
}
#===================================================================================================

1;
