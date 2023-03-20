use strict;
use warnings;

# as compared to other svX extractions, svLong:
#   expects PAF format as input
#   expects single-end long reads, so gaps are irrelevant
#   doesn't process outer clips, they are largely irrelevant and untrustworthy in likely being adapters and/or low quality bases
#   tracks strands, not sides, to allow path tracking across multiple SV junctions in a single molecule
#   has no sense of canonical strands, since an SV might occur in different parts/orientations of a path across multiple molecules
#   reports insertions as a positive _insertion size_, microhomology as negative (as opposed to inserttion = negative _overlap size_ elsewhere)

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
use vars qw($MIN_SV_SIZE @alnNodes @alnTypes @alnMapQs @alnSizes @alnInsSizes @alnAlns);   
my $minCigarSvDigits = length($MIN_SV_SIZE);

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
sub processAlignedSegment { 
    my ($aln) = @_;
    my ($cigar) = ($$aln[PAF_TAGS] =~ m/cg:Z:(\S+)/);

    # split a contiguous alignment conjoined around a smaller SV encoded in CIGAR by minimap2
    if($cigar =~ m/\d{$minCigarSvDigits}(D|I)/){        
        my $indexPos = $$aln[RSTART]; # RSTART is 0-based, on top genome strand, i.e., works L to R regardless of alignment strand
        my ($prevSize, $prevOp, $pendingOp) = (0, NULL_OP, NULL_OP);

        # step through all CIGAR operations; first and last operations are always M
        while ($cigar =~ (m/(\d+)(\w)/g)) { 
            my ($size, $operation) = ($1, $2);

            # handle largeD->smallI and largeD->largeI operations
            if($operation eq INSERTION and $pendingOp eq DELETION){
                $alnInsSizes[$#alnInsSizes] = $size;
                $pendingOp = NULL_OP;
            
            # handle largeI->smallD and largeI->largeD operations
            } elsif($operation eq DELETION and $pendingOp eq INSERTION){
                $alnTypes[$#alnTypes] = DELETION;
                $alnSizes[$#alnSizes] = $size;
                $indexPos += $size;
                $$aln[RSTART] = $indexPos; 
                $pendingOp = NULL_OP;

            # process a large indel not previoulsy handled in a sequential operation               
            } elsif($size >= $MIN_SV_SIZE and $operation ne MATCH_OP){
                my $isSmallDLargeI = ($prevOp eq DELETION and $operation eq INSERTION); # handle smallD->largeI operations

                # commit the alignment upstream of this large indel
                my @partial = @$aln; 
                $partial[REND] =  $isSmallDLargeI ? $indexPos - $prevSize : $indexPos;
                commitAlignmentNodes(\@partial);

                # commit the required nodes and junctions for this large indel
                if($operation eq INSERTION){
                    push @alnTypes,    $operation; 
                    push @alnMapQs,    0;
                    push @alnSizes,    $isSmallDLargeI ? $prevSize : 0;
                    push @alnInsSizes, $size; 
                    push @alnAlns,     [];
                } else {
                    $indexPos += $size;
                    push @alnTypes,    $operation; 
                    push @alnMapQs,    0;
                    push @alnSizes,    $size;
                    push @alnInsSizes, $prevOp eq INSERTION ? $prevSize : 0; # handle smallI->largeD operations
                    push @alnAlns,     [];
                }
                $pendingOp = $operation; # all paths except this one reset pendingOp to NULL_OP
                $$aln[RSTART] = $indexPos; 

            # jump over M and small D operations; isolated small insertions are ignored as they don't change indexPos
            } elsif($operation ne INSERTION) {
                $pendingOp = NULL_OP;
                $indexPos += $size;
            } else {
                $pendingOp = NULL_OP;
            }
            ($prevSize, $prevOp) = ($size, $operation);
        }
        commitAlignmentNodes($aln);

    # a single alignment with at most small indels not called as SVs
    } else {
        commitAlignmentNodes($aln);
    }

    # maintain proper 5'-3' node order on bottom strand
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
    push @alnInsSizes, 0;
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
        insSize => $$aln2[QSTART] - $$aln1[QEND] # i.e., microhomology is a negative number for svLong
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

1;
