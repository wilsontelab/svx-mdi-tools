use strict;
use warnings;

# STEP 2 - merge the two reads of a pair when they are overlapping, if not already merged by fastp

# constants (must match other scripts since not using perl package)
use constant {
    READ_PAIR_ID => 0, # columns in read-pair lines
    MOL_STRAND => 1,
    SEQ1 => 2, # SEQ indices refer to reads
    QUAL1 => 3,
    SEQ2 => 4,
    QUAL2 => 5,
    QNAME => 6,
    SEQ_MERGED => 7,
    QUAL_MERGED => 8,
    #-----------------
    MOL_MARKER => 0, # columns in molecule lines
    MOL_ID => 1,
    STRAND_COUNT1 => 2, # indices refer to strands
    STRAND_COUNT2 => 3,
    STRAND1 => 4, # indices refer to reads   
    STRAND2 => 5,
    UMI1 => 6,
    UMI2 => 7,
    IS_DUPLEX => 8,
    IS_MERGED => 9,
    #-------------------
    STRAND1 => 0, # consensus array indices, for code readability
    STRAND2 => 1, # so can have combinations like STRAND1|READ1|SEQ or DUPLEX|MERGED|QUAL
    DUPLEX => 2,    
    READ1 => 0,
    READ2 => 1,
    MERGED => 2,
    SEQ => 0,
    QUAL => 1,
};

# operating parameters
my $readLen = $ENV{READ_LEN};
my @mergeRegisters = (0,-1,1,-2,2,-3,3,-4,4); # ordered process of trying to merge proper pair overlaps
my $minMergeOverlapP  = $ENV{MIN_OVERLAP_PROPER} || 5;# the number of bases that must overlap between proper read ends for them to be merged
my $minMergeOverlapSV = $ENV{MIN_OVERLAP_SV} || 10;   # use a more stringent overlap requirement for SVs since not guided by alignment
our $minMergeDensity  = $ENV{MIN_MERGE_DENSITY} || 11/12; # at least this fraction of overlap bases must be matches to accept a merge
my $minMergeQual = $ENV{MIN_MERGE_QUAL} || 0.8; # reject merges of SV pairs with less than this MERGE_QUAL (max possible = 1)
my $bestMergeQual = 1; # MERGE_QUAL value for proper pairs, which were guided by genome alignment or fastp
our @maxMergeMismatch = map { $_ - int($_ * $minMergeDensity) - 1 } 0..($readLen*3);

# regular expressions
my $leftClip_  = qr/^(\d+)S/;
my $rightClip_ = qr/(\d+)S$/;

# working variables
use vars qw(@readPairs $mol @strands @downsample %ACGTNMatches %noIndelMatches);
my $NA = 0;
my $nullSymbol = '*';
my ($proxSeq, $distSeq, $proxQual, $distQual);

# execute merge action over all read pairs from a single source molecule
# MERGE_QUAL indicates confidence in the specificity of the merge
# sets $$mol[IS_MERGED] based on $minMergeDensity and $minMergeQual
sub mergeReads {
    
    # initialize the merging process
    @$mol[IS_MERGED, MERGE_QUAL] = (0, $NA);
    my ($mergeSub, @mergeParameters);
    
    # determine the starting process used to merge reads
    if($$mol[PAIR_OFFSET] eq $nullSymbol){
        $mergeSub = \&mergeWithoutGuidance; # TODO: is this even worth trying? fastp has already failed at merging...
        ($proxSeq, $distSeq, $proxQual, $distQual) = (SEQ1, SEQ2, QUAL1, QUAL2);        
    } else {        
        $mergeSub = \&mergeWithGuidance;
        if($$mol[PAIR_OFFSET] >= 0){ # reads potentially overalapping in their 3' ends
            ($proxSeq, $distSeq, $proxQual, $distQual) = (SEQ1, SEQ2, QUAL1, QUAL2);
            @mergeParameters = ($$mol[PAIR_OFFSET]); 
        } else { # reads sequenced past genomic into UMI+adapter (on left side at least)
            ($proxSeq, $distSeq, $proxQual, $distQual) = (SEQ2, SEQ1, QUAL2, QUAL1);
            @mergeParameters = (-$$mol[PAIR_OFFSET], 1);
        } 
    }
   
    # run merge process on all read pairs from a source molecule
    my ($refReadPair, %mergedPairs, %preMergedStrands);
    foreach my $strand(@strands){
        foreach my $readPair(@{$readPairs[$strand]}[@{$downsample[$strand]}]){
            
            # readPair merged previously by fastp, pass as is
            # merge parameters are already set based on this merge, if needed for other molecules in group
            if($$readPair[RP_IS_MERGED]){ 
                $$readPair[SEQ_MERGED]  = $$readPair[SEQ1];
                $$readPair[QUAL_MERGED] = $$readPair[QUAL1];
                $$readPair[MERGE_QUAL_RP] = $bestMergeQual;                
                push @{$mergedPairs{$strand}}, $readPair;
                !$refReadPair and $refReadPair = $readPair;
                $preMergedStrands{$strand}++;
                
            # otherwise, mergeSub sets values in readPair and returns guiding info for next pair    
            } else {
                my ($pairOffset, $isOverrunLeft, $mergeQual) = &$mergeSub($readPair, @mergeParameters);
                if(defined $pairOffset and $mergeQual >= $minMergeQual){ # seeding merge must be high quality
                    $mergeSub = \&mergeWithGuidance; # later pairs always have guidance from prior merges                   
                    !defined($mergeParameters[0]) and @mergeParameters = ($pairOffset, $isOverrunLeft);
                    push @{$mergedPairs{$strand}}, $readPair;
                    !$refReadPair and $refReadPair = $readPair;  
                }                  
            } 
        }
    }
    
    # coordinate merging between the two duplex strands (must be the same)
    my $nStrandsMerged = scalar keys %mergedPairs;
    if($nStrandsMerged == 0){
        return; # definitively unmerged on all paths to date
    } elsif($nStrandsMerged != @strands){ # two duplex strands are in conflict re: merging
        if(scalar keys %preMergedStrands == 0){
            return; # OK to continue with default unmerged state as all molecules came to us as pairs
        } else { # resolve this critical conflict by dropping the unmerged (likely low quality) strand
            @strands = keys %preMergedStrands;
        }
    }

    # report quality and length from the first encountered high quality merge
    # this is the pair that provided the guidance to other pairs
    $$mol[MERGE_QUAL] = @$refReadPair[MERGE_QUAL_RP];
    $$mol[IS_MERGED] = 1;
    
    # adjust downsample indices to reject any unmerged read pairs when some did merge
    foreach my $strand(@strands){
        @{$mergedPairs{$strand}} == @{$downsample[$strand]} and next; # all read pairs were merged
        my @is;
        foreach my $i(@{$downsample[$strand]}){
            my $rp = $readPairs[$strand][$i];
            $$rp[MERGE_QUAL_RP] and $$rp[MERGE_QUAL_RP] >= $minMergeQual and push @is, $i;
        }
        @{$downsample[$strand]} = @is;
    }
}

# merge two reads when we have guidance regarding the expected overlap register
# TODO: with improved guidance and fast pre-merging, this can become even more aggressive
sub mergeWithGuidance { 
    my ($readPair, $pairOffset, $isOverrunLeft) = @_;
    # ------------------     ----------------   proximal
    #    ------------------  -------------      distal; molDeltaLeft always >= 0

    # immediately reject if reads are not overlapping
    my $proxReadLen = length($$readPair[$proxSeq]);    
    !$isOverrunLeft and $pairOffset > $proxReadLen - $minMergeOverlapP and return (undef);
    my $distReadLen = length($$readPair[$distSeq]);    
    
    # try in predicted register first, but if fail also try shifted registers
    # stop trying once we have a good match, since we have guidance that it is the correct register    
    foreach my $register(@mergeRegisters){
        
        # get the candidate overlap based on the alignment positions
        my $wrkMolDeltaLeft = $pairOffset + $register;
        my $overlapLen = $proxReadLen - $wrkMolDeltaLeft;
        $overlapLen >= $minMergeOverlapP or next;
        $overlapLen <= $proxReadLen or next; # in case somebody else trimmed our reads, etc.
        $overlapLen <= $distReadLen or next;
        my $proxOverlap = substr($$readPair[$proxSeq], -$overlapLen);
        my $distOverlap = substr($$readPair[$distSeq], 0, $overlapLen);

        # check for overlap identity on the two reads
        if($proxOverlap eq $distOverlap){
            return stitchMerged($readPair, $isOverrunLeft, $wrkMolDeltaLeft, $proxOverlap, $overlapLen, $bestMergeQual); # mergeQual only used on first read pair              
        
        # if not identical, use slower process to match with attention to N bases
        } else {
            # TODO: here is where we could use e.g. Smith-Waterman to merge more aggressively
            my ($nMatch, $overlapSeq) = noIndelMatch($proxOverlap, $distOverlap);
            if($overlapSeq and $nMatch / $overlapLen >= $minMergeDensity){
                return stitchMerged($readPair, $isOverrunLeft, $wrkMolDeltaLeft, $overlapSeq, $overlapLen, $bestMergeQual);   
            } 
        }
    }
    return (undef);
}

# merge two reads WITHOUT the benefit of mapped positions or other guidance
# applies to the first merge attempt on anomalous read pairs
sub mergeWithoutGuidance { 
    my ($readPair, $isOverrunLeft) = @_;
    
    # score merges over all possible SV merge registers
    my %merges;
    my $proxReadLen = length($$readPair[$proxSeq]);
    my @oLens = $minMergeOverlapSV..min($proxReadLen, length($$readPair[$distSeq]));
    foreach my $overlapLen($isOverrunLeft ? reverse(@oLens) : @oLens){
        
        # get the candidate overlap based on the alignment positions
        my $proxOverlap = substr($$readPair[$proxSeq], -$overlapLen);
        my $distOverlap = substr($$readPair[$distSeq], 0, $overlapLen);
        
        # check for overlap identity on the two reads
        if($proxOverlap eq $distOverlap){
            $merges{$overlapLen}{$overlapLen} = $proxOverlap;              
        
        # if not identical, use slower process to match with attention to N bases
        } else {
            my ($nMatch, $overlapSeq) = noIndelMatch($proxOverlap, $distOverlap);
            if($overlapSeq and $nMatch / $overlapLen >= $minMergeDensity){
                $merges{$nMatch}{$overlapLen} = $overlapSeq;           
            } 
        }
    }
    
    # pick the best merge registers = shortest overlap from among the highest scoring overlaps
    # set the difference between the best and next overlaps as the MERGE_QUAL
    if(%merges){
        my @nMatches = sort {$b <=> $a} keys %merges;
        my @overlapLens = sort {$a <=> $b} keys %{$merges{$nMatches[0]}};
        my $mergeQual;
        if(@overlapLens > 1){
            $mergeQual = 0; # more than one merge register has the same nMatches
        } else {
            $mergeQual = $nMatches[0] - ($nMatches[1] || 0);
        }
        $mergeQual >= $minMergeQual and return stitchMerged( # only return high quality unguided merges
            $readPair, $isOverrunLeft,
            $proxReadLen - $overlapLens[0],
            $merges{$nMatches[0]}{$overlapLens[0]},
            $overlapLens[0],
            $mergeQual / $overlapLens[0]
        );
    }
    
    # if failed on 1st pass, retry assuming overrun into adapters
    if($isOverrunLeft){
        ($proxSeq, $distSeq, $proxQual, $distQual) = (SEQ1, SEQ2, QUAL1, QUAL2); # prep for next unguided read pair
        return (undef); # tried in all registers, report failure
    } else {
        ($proxSeq, $distSeq, $proxQual, $distQual) = (SEQ2, SEQ1, QUAL2, QUAL1);
        return mergeWithoutGuidance($readPair, 1);
    }
}

# assemble the final merged read sequences based on overlap sequence and register
# QUAL at this stage ~ probability that the called base represents the molecule seen by the sequencer
sub stitchMerged {
    my ($readPair, $isOverrunLeft, $pairOffset, $overlapSeq, $overlapLen, $mergeQual) = @_;
    if($isOverrunLeft){ # read pairs not also overrunning on right will fail merging and not get here
        $$readPair[SEQ_MERGED]  = $overlapSeq;
        $$readPair[QUAL_MERGED] = getConsensusQual_37($overlapSeq);
    } else {
        my $proxStemLen = length($$readPair[$proxSeq]) - $overlapLen;
        my $distStemLen = length($$readPair[$distSeq]) - $overlapLen;
        $$readPair[SEQ_MERGED]  = ($proxStemLen > 0 ? substr($$readPair[$proxSeq], 0, $proxStemLen) : "")
                                  .$overlapSeq.
                                  ($distStemLen > 0 ? substr($$readPair[$distSeq], -$distStemLen)   : "");
        $$readPair[QUAL_MERGED] = ($proxStemLen > 0 ? substr($$readPair[$proxQual], 0, $proxStemLen) : "") # pass non-confirmed tail bases as original qual
                                  .getConsensusQual_37($overlapSeq).
                                  ($distStemLen > 0 ? substr($$readPair[$distQual], -$distStemLen)   : "");
    }
    $$readPair[MERGE_QUAL_RP] = $mergeQual;
    return ($pairOffset, $isOverrunLeft, $mergeQual, $overlapLen);
}

1;

