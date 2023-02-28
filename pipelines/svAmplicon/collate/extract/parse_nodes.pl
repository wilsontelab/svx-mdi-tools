use strict;
use warnings;

# as compared to other svX extractions, svAmplicon
#   has very little need for processing outer clips; those positions are fixed as the amplicon ends
#   unmerged pair GAP junctions are assessed relative to the two outermost alignments to yield M=merge failue or N=inserted sequence

# constants
use constant {
    AMP_AMPLICON_ID => 0, # amplicon fields
    AMP_PROPER => 1,
    AMP_MOL_COUNT => 2,
    AMP_CHROM1 => 3,
    AMP_SIDE1 => 4,
    AMP_POS1 => 5,
    AMP_REF1 => 6,
    AMP_CHROM2 => 7,
    AMP_SIDE2 => 8,
    AMP_POS2 => 9,
    AMP_REF2 => 10,
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
    QUAL => 10,
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
    IS_REFERENCE => 3,
    MOL_COUNT => 4, 
    IS_MERGED => 5, 
    READ_N => 6, 
    MOL_CLASS => 6, # values added by extract_nodes regardless of pipeline or bam source
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
use vars qw($READ_LEN $MAX_INSERT_SIZE $MIN_SV_SIZE
            @alns @mol $amplicon);    

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
# source molecules that had one contiguous alignment, only possible in expectOverlap mode
sub commitContiguousAlignment { # aln1/umi1 are the left/proximal end of the source molecule as provided to the aligner
    my ($aln) = @_;
    $mol[MOL_CLASS] = IS_PROPER; # molecule may have an indel that did not cause alignment splitting, set in summarizeMolecule
    summarizeMolecule([$aln], "*"); # molecule has no inner endpoints
}
#---------------------------------------------------------------------------------------------------
# a merged read-pair (or single read) aligned in at least two segments, i.e., that crossed an SV junction
sub parseMergedSplit {
    # shortest clip1 will match umi1, longest will match umi2
    my ($outData1s, @alnIs) = sortReadAlignments(\@alns, LEFT, RIGHT);
    $mol[MOL_CLASS] = IS_SV;
    my $jxns = printSequencedJunctions(@alns[@alnIs]); 
    summarizeMolecule([@alns[@alnIs]], $jxns); # molecule has no inner endpoints
}
#---------------------------------------------------------------------------------------------------
# an unmerged FR read pair with no split alignments (i.e., exactly one aln each read)
# either a merge failure or anomalous due to an unaligned junction (possibly at an inner clip)
sub parseUnmergedHiddenJunction {
    my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
    $mol[MOL_CLASS] = IS_SV;
    setEndpointData(\my @innData1, $alns[$read1], RIGHT, LEFT);
    setEndpointData(\my @innData2, $alns[$read2], LEFT,  RIGHT); 
    my $jxnType = getJxnType($alns[$read1], $alns[$read2], \@innData1, \@innData2, GAP);
    my $jxn =  "*";
    if($jxnType eq PROPER or $jxnType eq MERGE_FAILURE){
        $mol[MOL_CLASS] = $jxnType; 
    } else {   
        $jxn = printJunction(GAP, $alns[$read1], $alns[$read2]);    
    }
    summarizeMolecule([$alns[$read1], $alns[$read2]], $jxn);  
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
            # @alns = @{$alnsByRead[READ1]}; # eliminate one read, it is identical to the other
            return;
        }  
    }
    # determine the order of the alignments along the sequenced molecule in the split reads
    foreach my $read(READ1, READ2){
        my ($fwdSide, $revSide) = ($read == READ1) ? (LEFT, RIGHT) : (RIGHT, LEFT);
        if(@{$alnsByRead[$read]} == 1){ # an unsplit read paired with a split read
            @{$alnIs[$read]} = (0); # the index within @alnsByRead, not @alns
        } else { # one of either one or two split reads
            ($outDatas[$read], @{$alnIs[$read]}) = sortReadAlignments($alnsByRead[$read], $fwdSide, $revSide); 
        }  
    }
    $mol[MOL_CLASS] = IS_SV;  
    # my $innI1 = $alnIs[READ1][0];
    # my $innI2 = $alnIs[READ2][0]; 
    my @jxns; 
    foreach my $read(READ1, READ2){
        push @jxns, @{$alnsByRead[$read]} > 1 ?
            printSequencedJunctions(@{$alnsByRead[$read]}[@{$alnIs[$read]}]) :
            "*";
        $read == READ1 and push @jxns, "GAP"; #printJunction(GAP, $alnsByRead[READ1][$innI1], $alnsByRead[READ2][$innI2])
    }    # as for all, aligns are ordered across the entire molecule, so alns[0] and alns[$#alns] are outemost
    summarizeMolecule(
        [@{$alnsByRead[READ1]}[@{$alnIs[READ1]}], @{$alnsByRead[READ2]}[reverse @{$alnIs[READ2]}]], 
        join(":::", @jxns)
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
    my @jxns;
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
        push @jxns, printJunction(SPLIT, $alns[READ1], $alns[READ2]);
    }
    join(":::", @jxns); # as throughout ":::" delimits elements order across the entire molecule
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
    
    # describe the junction
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $nodeClass);  
    my $svSize = getSvSize($innData[READ1], $innData[READ2], $jxnType);
    # my $isCanonical = isCanonicalStrand($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $jxnType); 

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

    # print rectified nodes and edges    
    my $overlap = "0/*"; # TODO: parse junction microhomology
    my $jxn = join(":", $nodeClass, $jxnType, $svSize, ${$innData[READ1]}[_NODE], ${$innData[READ2]}[_NODE], $overlap);
    isReportableSv($jxnType, $svSize) and printToJunctionFile($jxn);
    return $jxn;
}
sub printToJunctionFile {
    my ($jxn) = @_;
    # print join("\t", $jxn, @mol[AMPLICON_ID, MOL_ID, MOL_COUNT]), "\n";
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
    if($nodeClass == GAP){ # based on amplicon type
        if($$amplicon[AMP_PROPER] eq "expectOverlap"){ # read pair should have merged, why didn't it?
            $dist <= 0 and return MERGE_FAILURE; # these should be mostly low quality reads
            $dist  > 0 and return INSERTION;     # but novel sequence in the span could also lead to failed merging
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
sub getSvSize { # always a positive integer, zero if NA
    my ($innData1, $innData2, $jxnType) = @_;
    ($jxnType eq TRANSLOCATION or
     $jxnType eq UNKNOWN or 
     $jxnType eq INSERTION or 
     $jxnType eq PROPER) and return 0;
    $jxnType eq INVERSION and return abs($$innData2[_POS] - $$innData1[_POS]);
    $jxnType eq DELETION and return $$innData2[_POS] - $$innData1[_POS] - 1;
    ($jxnType eq DUPLICATION or
     $jxnType eq MERGE_FAILURE) and return $$innData1[_POS] - $$innData2[_POS] + 1;
    return 0;
}
sub isReportableSv {
    my ($jxnType, $svSize) = @_;
    ($jxnType eq TRANSLOCATION or
     $jxnType eq INSERTION) and return 1;
    ($jxnType eq INVERSION or
     $jxnType eq DELETION or 
     $jxnType eq DUPLICATION) and return $svSize >= $MIN_SV_SIZE;
     return 0;
}
# determine if an alignment pair is from the canonical strand of a junction fragment
# this is a junction-level property and thus could differ from the moleculeStrand
#   the canonical strand is:
#       the top genome strand for PDL and tranlocations handled as PDL
#       the top genome strand for the alignment NOT in the inverted segment
#           e.g., the top strand to the left of the left inversion junction and vice versa
# sub isCanonicalStrand {
#     my ($aln1, $aln2, $innData1, $innData2, $jxnType) = @_;
#     if($$innData1[_SIDE] eq $$innData2[_SIDE]){ # inversion-type handling, sort by pos along conjoined chromosomes
#         if($jxnType eq INVERSION){ # inversions
#             $$innData1[_POS] < $$innData2[_POS];
#         } else { # inversion-type translocations
#             $$aln1[RNAME_INDEX] < $$aln2[RNAME_INDEX];
#         }
#     } else { # del/dup and PDL-type translocations, canonical is LR pairs
#         $$innData1[_SIDE] eq LEFTWARD;
#     }
# }
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
    $$data[_NODE] = join("/", $$aln[RNAME], @$data[_SIDE, _POS]); # complete signature of the SV breakpoint position
}
#===================================================================================================

#===================================================================================================
# print molecule-defining data bits
#---------------------------------------------------------------------------------------------------
sub summarizeMolecule {
    my ($orderedAlns, $orderedJxns) = @_; 

    # calculate fragment size? (easy for merged, less clear for unmerged)
    # map map of reference

    # number of matched reference bases
    my $nMBases = 0; # as currently constructed, this includes segments inserted from outside the amplicon reference region
    map { while ($$_[CIGAR] =~ (m/(\d+)M/g)) { $nMBases += $1 } } @$orderedAlns;
    my $nDIBases = 0; # as currently constructed, this includes segments inserted from outside the amplicon reference region
    map { while ($$_[CIGAR] =~ (m/(\d+)[D|I]/g)) { $nDIBases += $1 } } @$orderedAlns;

    # calculate max internal MAPQ (requires 3 alignments at least)
    my $nAlns = scalar(@$orderedAlns);
    my $maxInternalMapq = "NA";
    if($nAlns >= 3){
        $maxInternalMapq = 0;
        map { $$_[MAPQ] > $maxInternalMapq and $maxInternalMapq = $$_[MAPQ] } @$orderedAlns[1..($nAlns-2)];
    }   

    # print one line per distinct molecule sequence
    print join("\t",
        $orderedJxns,
        join(":::", printMolCigars($orderedAlns)),
        join(":::", map { join(":", @$_[FLAG, RNAME, POS, MAPQ, CIGAR]) } @$orderedAlns),
        @mol, # molId ampliconId nOverlapBases isReference molCount merged MOL_CLASS 
        $nAlns,  
        $nMBases,
        $nDIBases,
        # $jxnTypes, # for filtering convenience     
        $maxInternalMapq,
        $$orderedAlns[0][SEQ],
        $$orderedAlns[0][QUAL],        
        $mol[IS_MERGED] ? "*" : $$orderedAlns[$nAlns - 1][SEQ],        
        $mol[IS_MERGED] ? "*" : $$orderedAlns[$nAlns - 1][QUAL]
    ), "\n"; 
}
#===================================================================================================

#===================================================================================================
# print large indels as trackable SVs
sub printMolCigars { # parse one or more alignment blocks per molecule, with large alignments joined by ":::"
    my ($orderedAlns) = @_;
    my @indels;
    foreach my $aln(@$orderedAlns){
        push @indels, printAlnCigar($$aln[RNAME], $$aln[POS], $$aln[CIGAR]);
    }
    @indels;
}
sub printAlnCigar { # one or more large CIGAR string indels in a single alignment block, joined by "::"
    my ($chrom, $start, $cigar) = @_;
    $cigar =~ s/\d+S//g;
    my $end = $start - 1;
    my @indels;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        if($1 >= $MIN_SV_SIZE and ($2 eq "I" or $2 eq "D")){
            my $indel = join(":",
                SPLIT, # i.e., the junction was sequenced outright
                $2 eq "I" ? "N" : "D", # can't use "I" for "insertion" since it is used for inversion SVs 
                $1, 
                join("/", $chrom, "L", $end),
                join("/", $chrom, "R", $end + $1 + 1),
                "?" # TODO: create process to examine for deletion microhomology, not immediately provided by aligner
            );
            printToJunctionFile($indel); # print each indel as a trackable junction for comparing molecules
            push @indels, $indel;
        }
        $2 eq "I" or $end += $1;
    }
    return @indels ? join("::", @indels): "*"; # here, "::" delimits multiple events within a single block (only possible for indels)
}
#===================================================================================================

1;
