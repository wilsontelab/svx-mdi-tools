use strict;
use warnings;

# as compared to other svX extractions, svAmplicon
#   has very little need for processing outer clips; those positions are fixed as the amplicon ends

# constants
use constant {
    AMP_AMPLICON_ID => 0, # amplicon fields
    AMP_PROPER => 1,
    AMP_MOL_COUNT => 2,
    AMP_CHROM1 => 3,
    AMP_SIDE1 => 4,
    AMP_POS1 => 5,
    AMP_REF1 => 6,
    AMP_PRIMER1 => 7,
    AMP_CHROM2 => 8,
    AMP_SIDE2 => 9,
    AMP_POS2 => 10,
    AMP_REF2 => 11,
    AMP_PRIMER2 => 12,
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
    HIDDEN_GAP => 0, # the gap junction when there are exactly two unmerged alignments
    SPLIT_GAP  => 1, # the gap junction when at least one unmerged read was split
    #-------------
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    PROPER        => "P",
    INSERTION     => "I", 
    MERGE_FAILURE => "M", # NOTE! these were changed for svAmplicon relative to other svx pipeline to account for small CIGAR events
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
    setEndpointData(\my @innData1, $alns[$read1], RIGHT, LEFT);
    setEndpointData(\my @innData2, $alns[$read2], LEFT,  RIGHT); 
    my $jxnType = getJxnType($alns[$read1], $alns[$read2], \@innData1, \@innData2, GAP, HIDDEN_GAP);
    my $jxn =  "*";
    if($jxnType eq PROPER or $jxnType eq MERGE_FAILURE){
        $mol[MOL_CLASS] = $jxnType; 
    } else {   
        $mol[MOL_CLASS] = IS_SV;
        $jxn = printJunction(GAP, HIDDEN_GAP, $alns[$read1], $alns[$read2]);    
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
        my ($fwdSide, $revSide) = ($read == READ1) ? (LEFT, RIGHT) : (RIGHT, LEFT); # thus, sorts outer to inner
        if(@{$alnsByRead[$read]} == 1){ # an unsplit read paired with a split read
            @{$alnIs[$read]} = (0); # the index within @alnsByRead, not @alns
        } else { # one of either one or two split reads
            ($outDatas[$read], @{$alnIs[$read]}) = sortReadAlignments($alnsByRead[$read], $fwdSide, $revSide); 
        }  
    }
    $mol[MOL_CLASS] = IS_SV;  
    my $innI1 = $alnIs[READ1][$#{$alnIs[READ1]}]; # the innermost alignments over each read
    my $innI2 = $alnIs[READ2][$#{$alnIs[READ2]}]; 
    my @jxns; 
    foreach my $read(READ1, READ2){
        push @jxns, @{$alnsByRead[$read]} > 1 ?
            printSequencedJunctions(@{$alnsByRead[$read]}[@{$alnIs[$read]}]) :
            "*";
        $read == READ1 and push @jxns, printJunction(GAP, SPLIT_GAP, $alnsByRead[READ1][$innI1], $alnsByRead[READ2][$innI2]);
    }    
    summarizeMolecule( # as for all, final aligns are ordered across the entire molecule, so alns[0] and alns[$#alns] are outemost
        [@{$alnsByRead[READ1]}[        @{$alnIs[READ1]}], 
         @{$alnsByRead[READ2]}[reverse @{$alnIs[READ2]}]], 
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
        push @jxns, printJunction(SPLIT, undef, $alns[READ1], $alns[READ2]);
    }
    join(":::", @jxns); # as throughout ":::" delimits elements order across the entire molecule
}
# function for printing junctions in proper order
# here, outer and inner clips are for the _alignments_ immediately flanking the predicted junction
sub printJunction { 
    my ($nodeClass, $gapType, @alns) = @_;

    # collect endpoint data on the alignments that flank the junction
    # this step assigns node alignment sides based on FF strands
    # from this point forward we use node sides, not FLAG REVERSE, to track orientation
    my (@innData);
    setEndpointData(\@{$innData[READ1]}, $alns[READ1], RIGHT, LEFT);
    setEndpointData(\@{$innData[READ2]}, $alns[READ2], LEFT,  RIGHT);
    
    # describe the junction
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2], $nodeClass, $gapType);  
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

    # print rectified junctions  
    my ($overlap, $jxnBases) = getJxnStructure($nodeClass, $alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2]); 
    join(":", $nodeClass, $jxnType, $svSize, ${$innData[READ1]}[_NODE], ${$innData[READ2]}[_NODE], $overlap, $jxnBases);
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
    my ($aln1, $aln2, $innData1, $innData2, $nodeClass, $gapType) = @_;  

    # simplest case, translocation = alignment to different chromosomes
    # PDL vs. inversion sub-handling determined by isCanonicalStrand
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;

    # inversions determined by unexpected alignment orientiations
    $$innData1[_SIDE] eq $$innData2[_SIDE] and return INVERSION;

    # PDL are distinguish based on their position along the conjoined chromosomes
    my $dist = $$innData1[_SIDE] eq LEFTWARD ? 
                    $$innData2[_POS] - $$innData1[_POS] : # LR pair, canonical PDL strand
                    $$innData1[_POS] - $$innData2[_POS];  # RL pair, non-canonical strand (will be flipped later)
    if($nodeClass == GAP and $gapType == HIDDEN_GAP){ # the outer alignments of a simple unmerged pair have a predictable relationship based on amplicon design
        if($$amplicon[AMP_PROPER] eq "expectOverlap"){ # read pair should have merged, why didn't it?
            $dist <= 0 and return MERGE_FAILURE; # these should be mostly low quality reads
            $dist  > 0 and return INSERTION;     # but novel sequence in the span could also lead to failed merging
        } elsif($$amplicon[AMP_PROPER] eq "expectGaps") { # POOR DESIGN CHOICE, ignore until someone needs this path
            return UNKNOWN;                                # most will be proper anyway
            # my $refSpan = $$amplicon[AMP_POS2] - $$amplicon[AMP_POS1] + 1; # this code is a start but not adequate
            # my $expectedDist = $refSpan - 2 * $READ_LEN + 1;
            # my $gapDelLimit = $MAX_INSERT_SIZE - 2 * $READ_LEN;
            # $dist < $expectedDist and return DELETION;
            # $dist > $expectedDist and return INSERTION;
        } else { # any outer alignment pairs from "notPossible" amplicons are SVs of the type consistent with the amplicon design
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
sub isReportableSv { # TODO: unused here, move to script after molecules types are grouped by junction sequence
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
sub getJxnStructure {
    my ($nodeClass, $aln1, $aln2, $innData1, $innData2) = @_;
    $nodeClass eq GAP and return ("NA", "NA");
    my $readPos1 = length($$aln1[SEQ]) - $$innData1[_CLIP];
    my $readPos2 = $$innData2[_CLIP] + 1;
    my $overlap = $readPos1 - $readPos2 + 1;
    my $jxnBases = "*"; # a blunt joint
    if($overlap > 0){  # microhomology
        $jxnBases = substr($$aln1[SEQ], $readPos2,     $overlap); # TODO: any value to reverse complementing inversions?
    } elsif($overlap < 0){ # inserted bases
        $jxnBases = substr($$aln1[SEQ], $readPos1 + 1, $overlap);
    }
    ($overlap, $jxnBases);
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
    $$data[_NODE] = join("/", $$aln[RNAME], @$data[_SIDE, _POS]); # complete signature of the SV breakpoint position
}
#===================================================================================================

#===================================================================================================
# print molecule-defining data bits
#---------------------------------------------------------------------------------------------------
sub summarizeMolecule {
    my ($orderedAlns, $orderedJxns) = @_; 

    # exact or estimated size of the sequenced DNA molecule
    my $nAlns = scalar(@$orderedAlns);    
    my $insertSize = $mol[IS_MERGED] ? 
        length($$orderedAlns[0][SEQ]) :
        ($mol[N_OVERLAP_BASES] eq "NA" ? 
            "NA" : 
            length($$orderedAlns[0][SEQ]) + length($$orderedAlns[$nAlns -1][SEQ]) - $mol[N_OVERLAP_BASES]
        );

    # number of M, D and I reference bases
    my $nMBases = 0; # as constructed, this includes segments inserted from outside the amplicon reference region
    map { while ($$_[CIGAR] =~ (m/(\d+)M/g)) { $nMBases += $1 } } @$orderedAlns;
    my $nDBases = 0; # as currently constructed, this includes segments inserted from outside the amplicon reference region
    map { while ($$_[CIGAR] =~ (m/(\d+)D/g)) { $nDBases += $1 } } @$orderedAlns;
    my $nIBases = 0; # as currently constructed, this includes segments inserted from outside the amplicon reference region
    map { while ($$_[CIGAR] =~ (m/(\d+)I/g)) { $nIBases += $1 } } @$orderedAlns;

    # max internal alignment MAPQ (requires 3 alignments at least)
    my $maxInternalMapq = "NA";
    if($nAlns >= 3){
        $maxInternalMapq = 0;
        map { $$_[MAPQ] > $maxInternalMapq and $maxInternalMapq = $$_[MAPQ] } @$orderedAlns[1..($nAlns-2)];
    }  

    # parse an abbreviated junction type sequence suitable for tabulating and filtering
    my @indels = printMolCigars($orderedAlns);
    my @orderedJxns = split(":::", $orderedJxns);
    my @jxnsShort;    
    foreach my $i(0..($nAlns - 1)){
        if($indels[$i] ne "*"){
            foreach my $indel(split("::", $indels[$i])){
                my ($nodeClass, $jxnType) = split(":", $indel, 3);
                $jxnType and push @jxnsShort, lc($jxnType); # CIGAR indels are reported as lower case (implies smaller)
            }
        }
        if($orderedJxns[$i]){
            my ($nodeClass, $jxnType) = split(":", $orderedJxns[$i], 3);
            $jxnType and push @jxnsShort, $jxnType;
        }
    } 

    # print one line per distinct molecule sequence, packed for efficient sorting and grouping
    print join("\t",

        # first three columns are used to group molecules by SV path through molecule
        $mol[AMPLICON_ID],
        $orderedJxns,
        join(":::", @indels), 

        # next two columns are used to select the best molecule to represent the SV path 
        @mol[IS_REFERENCE, MOL_COUNT],

        # remaining columns describe the molecule in detail
        # as molecules are grouped, columns may be aggregated as first, mean, etc.
        join(":::", map { join(":", @$_[FLAG, RNAME, POS, MAPQ, CIGAR]) } @$orderedAlns),
        @mol[MOL_ID, MOL_CLASS, IS_MERGED, N_OVERLAP_BASES],
        $insertSize, # describe the span of molecule bases
        $nMBases,
        $nDBases,
        $nIBases,
        int(
            ($mol[IS_MERGED] ? 
                getAvgQual($$orderedAlns[0][QUAL]) : 
                getAvgQual($$orderedAlns[0][QUAL].$$orderedAlns[$nAlns - 1][QUAL])
            ) + 0.5
        ),
        $nAlns, # describe the span of molecule alignments and junctions
        scalar(@jxnsShort),
        @jxnsShort ? join("", @jxnsShort) : "*",
        $maxInternalMapq,
        $$orderedAlns[0][SEQ], # provide the SEQ and QUAL of the index molecule
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
        if($1 >= $MIN_SV_SIZE and ($2 eq INSERTION or $2 eq DELETION)){
            push @indels, join(":",
                SPLIT, # i.e., the junction was sequenced outright
                $2,
                $1, 
                join("/", $chrom, "L", $end),
                join("/", $chrom, "R", $end + $1 + 1),
                "NA", # TODO: create process to examine for deletion microhomology, not immediately provided by aligner
                "NA"
            );
        }
        $2 eq INSERTION or $end += $1;
    }
    return @indels ? join("::", @indels): "*"; # here, "::" delimits multiple events within a single block (only possible for indels)
}
#===================================================================================================

1;
