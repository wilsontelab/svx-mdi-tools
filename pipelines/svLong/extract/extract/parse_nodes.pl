use strict;
use warnings;

# as compared to other svX extractions, svLong
#   expects PAF format as input
#   expect single-ended long reads, so gaps are irrelevant
#   doesn't process outer clips, they are largely insignificant and untrustworthy in likely being low quality bases
#   doesn't provide sequence-level analysis, since ONT quality is not sufficient on its own

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
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    ALN1 => 0,
    ALN2 => 1,
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
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    CIGAR_MATCH   => "M",
    #-------------
    _NULL => '*',   
};

# # operating parameters
# sub getLeftClip  { $_[0][QSTART] } 
# sub getRightClip { $_[0][QLEN] - $_[0][QEND] }
# my @clips = (\&getRightClip, \&getLeftClip);
# my @sides = (RIGHTWARD, LEFTWARD); # e.g., at a left end clip, the read continues righward from that point

# working variables
use vars qw($MIN_SV_SIZE $WINDOW_POWER $WINDOW_SIZE 
            @alnNodes @alnTypes @alnMapQs @alnSizes @alnInsSizes 
            $molId);    
my $minCigarSvDigits = $WINDOW_POWER + 1; # so, 10**3 == 1000 requires 4-digit CIGAR operation to be called an SV
my $minSvSize = $WINDOW_SIZE + 1; # thus, indel operation is guaranteed to cross a window border

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
sub processAlignedSegment { 
    my ($aln) = @_;
    my ($cigar) = ($$aln[PAF_TAGS] =~ m/cg:Z:(\S+)/);

    # split a contiguous alignment conjoined around a smaller SV encoded in CIGAR by minimap2
    if($cigar =~ m/\d{$minCigarSvDigits}(D|I)/){        
        my $indexPos = $$aln[RSTART]; # RSTART is 0-based, on top genome strand, i.e., works L to R regardless of alignment strand
        my ($prevSize, $prevOp, $pendingOp) = (0, "X", "X");

        # step through all CIGAR operations; first and last operations are always M
        while ($cigar =~ (m/(\d+)(\w)/g)) { 
            my ($size, $operation) = ($1, $2);

            # handle largeD->smallI and largeD->largeI operations
            if($operation eq INSERTION and $pendingOp eq DELETION){
                $alnInsSizes[$#alnInsSizes] = $size;
                $pendingOp = "X";
                print join(" ", $prevSize, $prevOp, $size, $operation, "largeD->anyI"), "\n";
            
            # handle largeI->smallD and largeI->largeD operations
            } elsif($operation eq DELETION and $pendingOp eq INSERTION){
                $alnTypes[$#alnTypes] = DELETION;
                $alnSizes[$#alnSizes] = $size;
                $indexPos += $size;
                $$aln[RSTART] = $indexPos; 
                $pendingOp = "X";
                print join(" ", $prevSize, $prevOp, $size, $operation, "largeI->anyD"), "\n";

            # process a large indel not previoulsy handled in a sequential operation               
            } elsif($size > $minSvSize and $operation ne CIGAR_MATCH){
                my $isSmallDLargeI = ($prevOp eq DELETION and $operation eq INSERTION); # handle smallD->largeI operations

                print join(" ", $prevSize, $prevOp, $size, $operation, $isSmallDLargeI ? "smallD->largeI" : ""), "\n";

                # commit the alignment upstream of this large indel
                my @partial = @$aln; 
                $partial[REND] =  $isSmallDLargeI ? $indexPos - $prevSize : $indexPos;
                incrementWindowCoverage(\@partial); 
                commitAlignmentNodes(\@partial);

                # commit the required nodes and junctions for this large indel
                if($operation eq INSERTION){
                    push @alnTypes,    $operation; 
                    push @alnMapQs,    0;
                    push @alnSizes,    $isSmallDLargeI ? $prevSize : 0;
                    push @alnInsSizes, $size; 

                } else {
                    $indexPos += $size;
                    push @alnTypes,    $operation; 
                    push @alnMapQs,    0;
                    push @alnSizes,    $size;
                    push @alnInsSizes, $prevOp eq INSERTION ? $prevSize : 0; # handle smallI->largeD operations
                }
                $pendingOp = $operation;
                $$aln[RSTART] = $indexPos; 

            # jump over M and small D operations; isolated small insertions are ignored
            } elsif($operation ne INSERTION) {
                $pendingOp = "X";
                $indexPos += $size;
            } else {
                $pendingOp = "X";
            }
            ($prevSize, $prevOp) = ($size, $operation);
        }
        incrementWindowCoverage($aln);
        commitAlignmentNodes($aln);

    # a single alignment with at most small indels
    } else {
        incrementWindowCoverage($aln);
        commitAlignmentNodes($aln);
    }

    # maintain proper 5'-3' node order on bottom strand
    $$aln[STRAND] eq "-" and @alnNodes = reverse @alnNodes;   
}
sub commitAlignmentNodes {
    my ($aln) = @_;
    push @alnNodes, (
        getSignedWindow(@$aln[RNAME_INDEX, RSTART, STRAND], 1),
        getSignedWindow(@$aln[RNAME_INDEX, REND,   STRAND], 0)
    );   
    push @alnTypes,    ALIGNMENT; 
    push @alnMapQs,    $$aln[MAPQ];
    push @alnSizes,    $$aln[REND] - $$aln[RSTART];
    push @alnInsSizes, 0;
}
#---------------------------------------------------------------------------------------------------
sub printSplitJunction { 
    my ($aln1, $aln2) = @_;
    my $nodePos1 = $$aln1[STRAND] eq "+" ? $$aln1[REND] : $$aln1[RSTART] + 1;
    my $nodePos2 = $$aln2[STRAND] eq "-" ? $$aln2[REND] : $$aln2[RSTART] + 1;
    my $jxnType = getJxnType($aln1, $aln2, $nodePos1, $nodePos2);  
    my $svSize = getSvSize($jxnType, $nodePos1, $nodePos2);    
    my $overlap = $$aln1[QEND] - $$aln2[QSTART] + 1;
    push @alnTypes,    $jxnType; 
    push @alnMapQs,    0;
    push @alnSizes,    $svSize;
    push @alnInsSizes, -$overlap; # i.e., microhomology is a negative number for svLong
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
    my ($aln1, $aln2, $nodePos1, $nodePos2) = @_;  
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;
    $$aln1[STRAND] ne $$aln2[STRAND] and return INVERSION;
    my $dist = $$aln1[STRAND] eq "+" ? 
        $nodePos2 - $nodePos1 : 
        $nodePos1 - $nodePos2;  
    $dist <= 0 and return DUPLICATION;
    $dist  > 1 and return DELETION; # pos delta is 1 for a continuous proper alignment
    return UNKNOWN; # should not occur, aligner should have made contiguoues, consider it an error condition
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
