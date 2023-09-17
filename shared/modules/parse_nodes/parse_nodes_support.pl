use strict;
use warnings;

# common functions to svAmplicon and svPore

# constants
use constant {
    DELETION      => "D",
    INSERTION     => "I", 
    MATCH_OP      => "M",
    NULL_OP       => "X"
};

# working variables
use vars qw($MIN_SV_SIZE $APPEND_JXN_BASES
            @alnNodes @alnMapQs @alnCigars @alnAlnQs @alnTypes @alnSizes @alnInsSizes @alnAlns
            @nodes    @mapQs    @cigars    @alnQs    @types    @sizes    @insSizes    @outAlns);  
my $minCigarSvDigits = length($MIN_SV_SIZE);
my %qryOps = map { $_ => 1 } qw(S H I M);

# step through a single alignment to find small SVs encoded in the CIGAR string
# called by relevant parse_nodes.pl scripts
sub parseSvsInCigar { 
    my ($aln, $refStartI, $refEndI, $cigar, $commitAlnFn, @extra) = @_;

    # step through each CIGAR operation of sufficient size
    if($cigar =~ m/\d{$minCigarSvDigits}(D|I)/){        
        my ($prevSize, $prevOp, $pendingOp, $alnCigar) = (0, NULL_OP, NULL_OP, "");

        # refPos is 0-based, on top genome strand, i.e., works L to R regardless of alignment strand
        my $refPos = $$aln[$refStartI];
        my $qryPos = 0;    

        # step through all CIGAR operations
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
                $refPos += $size;
                $$aln[$refStartI] = $refPos; 
                $pendingOp = NULL_OP;

            # process a large indel not previously handled in a sequential operation   
            } elsif($size >= $MIN_SV_SIZE and ($operation eq DELETION or $operation eq INSERTION)){

                # handle smallD->largeI operations
                my $isSmallDLargeI = ($prevOp eq DELETION and $operation eq INSERTION); 

                # commit any alignment upstream of this large indel
                my @partial = @$aln; 
                $partial[$refEndI] = $isSmallDLargeI ? $refPos - $prevSize : $refPos;
                &$commitAlnFn(\@partial, $alnCigar, @extra);

                # commit the required edges for this large indel
                if($operation eq INSERTION){
                    push @alnMapQs,    0;
                    push @alnCigars,   "NA";
                    push @alnAlnQs,    0;
                    push @alnTypes,    $operation; 
                    push @alnSizes,    $isSmallDLargeI ? $prevSize : 0;
                    push @alnInsSizes, $size.($APPEND_JXN_BASES ? "\tNA:$qryPos" : ""); # svAmplicon at least expect jxnBases appended to insSize; NA jxnBases flags CIGAR indels
                    push @alnAlns,     [];
                } else {
                    $refPos += $size;
                    push @alnMapQs,    0;
                    push @alnCigars,   "NA";
                    push @alnAlnQs,    0;
                    push @alnTypes,    $operation; 
                    push @alnSizes,    $size;
                    push @alnInsSizes, ($prevOp eq INSERTION ? $prevSize : 0).($APPEND_JXN_BASES ? "\tNA:$qryPos" : ""); # handle smallI->largeD operations
                    push @alnAlns,     [];
                }

                # all paths except this one reset pendingOp to NULL_OP
                $pendingOp = $operation; 

                # reset to the next split alignment
                $$aln[$refStartI] = $refPos; 
                $alnCigar = "";                 

            # jump over M and small D operations; 
            } elsif($operation eq MATCH_OP or $operation eq DELETION) {
                $refPos += $size;                
                $pendingOp = NULL_OP;
                $alnCigar .= "$size$operation";

            # clips and isolated small insertions are ignored as they don't change refPos    
            } else {
                $pendingOp = NULL_OP;
                $alnCigar .= "$size$operation";
            }
            ($prevSize, $prevOp) = ($size, $operation);
            $qryOps{$operation} and $qryPos += $size;
        }
        $alnCigar and &$commitAlnFn($aln, $alnCigar, @extra);

    # a single alignment with at most small indels not called as SVs
    } else {
        &$commitAlnFn($aln, $cigar, @extra);
    }
}

# set junction MAPQ as minimum MAPQ of the two flanking alignments
sub fillJxnQs {
    if(@mapQs > 1){
        for (my $i = 1; $i <= $#mapQs - 1; $i += 2){
            $mapQs[$i] = min($mapQs[$i - 1], $mapQs[$i + 1]);
            $alnQs[$i] = min($alnQs[$i - 1], $alnQs[$i + 1]);
        }
    }
}

# print a validated molecule to STDOUT
sub printMolecule {
    my ($molId, @extraVals) = @_;
    my $out = ""; # "----------------------------\n";  
    foreach my $i(0..$#types){
        $out .= join(
            "\t", 
            $molId,
            $nodes[$i], 
            $nodes[$i + 1], 
            $mapQs[$i], 
            $cigars[$i], 
            $alnQs[$i], 
            $types[$i], 
            $sizes[$i], 
            $insSizes[$i],
            @extraVals # nStrands for svPore, various for svAmplicon
        )."\n";
    }
    print $out; 
}

1;
