use strict;
use warnings;

# as compared to other svX extractions, svPore:
#   expects PAF format as input
#   expects single-end long reads, so gaps are irrelevant
#   doesn't process outer clips, they are largely irrelevant and untrustworthy in being adapters and/or low quality bases
#   tracks strands, not sides, to allow path tracking across multiple SV junctions in a single molecule
#   has no sense of canonical strands, since an SV might occur in different parts/orientations of a path across multiple molecules
#   reports insertions as a positive _insertion size_, microhomology as negative (as opposed to insertion = negative _overlap size_ elsewhere)

# constants
use constant {
    QNAME => 0, # PAF fields
    QLEN => 1,
    QSTART => 2,
    QEND => 3,
    STRAND => 4,
    RNAME => 5,
    RLEN => 6,
    RSTART => 7,
    REND => 8,
    N_MATCHES => 9,
    N_BASES => 10,
    MAPQ => 11,
    PAF_TAGS => 12,
    RNAME_INDEX => 13,  # added by us  
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    #-------------
    MATCH_OP      => "M",
    NULL_OP       => "X"
};

# working variables
use vars qw($MIN_SV_SIZE 
            @alnNodes @alnTypes @alnMapQs @alnSizes @alnInsSizes @alnAlns
            @nodes @types @mapQs @sizes @insSizes @outAlns);   
my $minCigarSvDigits = length($MIN_SV_SIZE);

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
sub processAlignedSegment { 
    my ($aln) = @_;

    # run the CIGAR string
    my ($cigar) = ($$aln[PAF_TAGS] =~ m/cg:Z:(\S+)/);
    parseSvsInCigar($aln, RSTART, REND, $cigar, \&commitAlignmentNodes);

    # maintain proper 5'-3' node order on bottom strand, since svPore tracks genome paths
    if($$aln[STRAND] eq "-"){
        @alnNodes    = reverse @alnNodes; 
        @alnTypes    = reverse @alnTypes;
        @alnMapQs    = reverse @alnMapQs;
        @alnSizes    = reverse @alnSizes;
        @alnInsSizes = reverse @alnInsSizes;
        @alnAlns     = reverse @alnAlns;
    }
}
sub commitAlignmentNodes { # add continguous "A" alignment segment to molecule chain
    my ($aln) = @_;        # often includes small indels in the alignment span
    push @alnNodes, (
        getSignedWindow(@$aln[RNAME_INDEX, RSTART, STRAND], 1),
        getSignedWindow(@$aln[RNAME_INDEX, REND,   STRAND], 0)
    );   
    push @alnTypes,    ALIGNMENT; 
    push @alnMapQs,    $$aln[MAPQ];
    push @alnSizes,    $$aln[REND] - $$aln[RSTART];
    push @alnInsSizes, join(
        "\t", 
        0,                         # insSize not applicable for alignments
        $$aln[QSTART], $$aln[QEND] # alignments carry query positions in xStart and xEnd
    );
    push @alnAlns,     $aln;
}
#---------------------------------------------------------------------------------------------------
sub processSplitJunction { 
    my ($aln1, $aln2) = @_;
    my $nodePos1 = $$aln1[STRAND] eq "+" ? $$aln1[REND] : $$aln1[RSTART] + 1;
    my $nodePos2 = $$aln2[STRAND] eq "-" ? $$aln2[REND] : $$aln2[RSTART] + 1;
    my $jxnType = getJxnType($aln1, $aln2, $nodePos1, $nodePos2);  
    {
        jxnType => $jxnType,
        svSize  => getSvSize($jxnType, $nodePos1, $nodePos2),
        insSize => join(
            "\t", 
            $$aln2[QSTART] - $$aln1[QEND],  # i.e., microhomology is a negative number for svPore
            $nodePos1, $nodePos2            # junctions carry reference positions in xStart and xEnd
        )
    }
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
sub getJxnType {
    my ($aln1, $aln2, $nodePos1, $nodePos2) = @_;  
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;
    $$aln1[STRAND] ne $$aln2[STRAND] and return INVERSION;
    my $dist = $$aln1[STRAND] eq "+" ? 
        $nodePos2 - $nodePos1 : 
        $nodePos1 - $nodePos2;  
    $dist <= 0 and return DUPLICATION;
    $dist  > 1 and return DELETION; # pos delta is 1 for a continuous proper alignment
    return UNKNOWN; # should not occur, aligner should have made a contiguous alignment, consider it an error condition
}
sub getSvSize { # always a positive integer, zero if NA
    my ($jxnType, $nodePos1, $nodePos2) = @_;
    ($jxnType eq TRANSLOCATION or
     $jxnType eq UNKNOWN or 
     $jxnType eq ALIGNMENT) and return 0;
    abs($nodePos2 - $nodePos1);
}
#===================================================================================================

#===================================================================================================
# check whether a molecule is consistent with a duplex foldback inversion
# if yes, keep first half only and mark molecules as duplex
# is most sensitive and acceptable to reject any inversion with reverse-complement overlap between its flanking alignments
#   ----->
#         | V
#   <-----
# strictly speaking we expect symmetry, but cannot count on complete alignment of both flanks
# we might expect adapters in the junction, but exploration says they may not be there
# and foldback inversions as so supsect as artifacts that we maintain high sensivitiy for purging them
#---------------------------------------------------------------------------------------------------
sub checkForDuplex {
    my ($nEdges) = @_;
    $nEdges > 1 or return 1;
    for(my $i = 0; $i <= $#types; $i++){
        $types[$i] eq INVERSION or next;     
        if($outAlns[$i - 1][RSTART] <= $outAlns[$i + 1][REND] and
           $outAlns[$i + 1][RSTART] <= $outAlns[$i - 1][REND]){
            my $j = $i - 1;
            @types = @types[0..$j]; # no need to splice @nodes, they'll never print
            @mapQs = @mapQs[0..$j];
            @sizes = @sizes[0..$j];
            @insSizes = @insSizes[0..$j];
            return 2; # i.e., two strands = duplex
        }
    }
    return 1; # i.e., one strand = simplex
}
#===================================================================================================

1;
