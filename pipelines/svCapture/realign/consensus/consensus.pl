use strict;
use warnings;

# STEP 3 - make a single-strand consensus for each read, or merged read, of each source strand
# STEP 4 - make a duplex consensus from the single-strand consensuses when available

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
    UMI1 => 4,
    UMI2 => 5,
    IS_DUPLEX => 6,
    IS_MERGED => 7,
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
my $minScoreFactor  = $ENV{MIN_SCORE_FACTOR} || 0.7; # e.g. 100 length * 0.7 factor = minimum acceptable _score_ of 70 to accept a consensus alignment
my $consensusFactor = $ENV{CONSENSUS_FACTOR} || 2/3; # the fraction of reads/bases that must have a common value to be called a consensus

# working variables
use vars qw(@readPairs @downsample @consensus %mOperations);
  
# make a consensus independently on the two unmerged reads of a pair
sub makeStrandReadConsensus {
    my ($strand) = @_;
    
    # source molecule+strand had only one unmerged read pair, use as is
    if(@{$readPairs[$strand]} == 1){
        $consensus[$strand][READ1][SEQ]  = ${$readPairs[$strand]}[0][SEQ1];
        $consensus[$strand][READ2][SEQ]  = ${$readPairs[$strand]}[0][SEQ2];
        $consensus[$strand][READ1][QUAL] = ${$readPairs[$strand]}[0][QUAL1]; # pass non-confirmed bases as original qual
        $consensus[$strand][READ2][QUAL] = ${$readPairs[$strand]}[0][QUAL2]; # don't downgrade quality yet, do that later  
    } else {

        # count all of the unique sequences for each read of the downsampled set of read-pairs
        my @readSeqs;
        foreach my $readPair(@{$readPairs[$strand]}[@{$downsample[$strand]}]){
            push @{$readSeqs[READ1]{$$readPair[SEQ1]}}, $$readPair[QUAL1] ;
            push @{$readSeqs[READ2]{$$readPair[SEQ2]}}, $$readPair[QUAL2] ; 
        }

        # determine the multi-sequence consensus per read
        foreach my $read(READ1, READ2){
            makeConsensusMulti($strand, $read, $readSeqs[$read]);
        }                  
    }
}

# make a strand consensus on the merged reads of a pair
sub makeStrandMergeConsensus {
    my ($strand) = @_;
    
    # source molecule+strand had only one merged read pair, use as is
    if(@{$readPairs[$strand]} == 1){
        $consensus[$strand][MERGED][SEQ]  = ${$readPairs[$strand]}[0][SEQ_MERGED]; # see comments above
        $consensus[$strand][MERGED][QUAL] = ${$readPairs[$strand]}[0][QUAL_MERGED];
    } else {
        
        # count all of the unique sequences for merged read of the downsampled set of read-pairs
        my %mergedSeqs;
        foreach my $readPair(@{$readPairs[$strand]}[@{$downsample[$strand]}]){
            push @{$mergedSeqs{$$readPair[SEQ_MERGED]}}, $$readPair[QUAL_MERGED] ;
        }
        
        # determine the multi-sequence consensus of the merged read
        makeConsensusMulti($strand, MERGED, \%mergedSeqs);   
    }
}

# make a strand sequence consensus from two or more input read pairs
# either on one read of the pair, or the merged reads
sub makeConsensusMulti {
    my ($strand, $read, $seqs) = @_;
    
    # source molecule+strand+read/merged had only one sequence value over multiple read-pairs, use as is
    my $nSeqVals = keys %$seqs;
    if($nSeqVals == 1){
        my $seq = (keys %$seqs)[0];
        $consensus[$strand][$read][SEQ]  = $seq;
        $consensus[$strand][$read][QUAL] = getConsensusQual_38($seq); 
        
    # exactly two reads that had different sequences; cannot resolve, mask the raw base ambiguities
    } elsif(@{$readPairs[$strand]} == 2){   
        $consensus[$strand][$read][SEQ]  = maskTwoSequences(keys %$seqs);
        $consensus[$strand][$read][QUAL] = getConsensusQual_38($consensus[$strand][$read][SEQ]);        
    } else {
        
        # sort multiple sequence values by their number of occurences
        my @sortedSeq = sort { @{$$seqs{$b}} <=> @{$$seqs{$a}} } keys %$seqs;        

        # source molecule+strand+read/merged had one predominant sequence (would win any base consensus counting)
        if(@{$$seqs{$sortedSeq[0]}} >= $nSeqVals * $consensusFactor){
            $consensus[$strand][$read][SEQ]  = $sortedSeq[0];    

        # complex sequence pattern requires base-by-base pattern evaluation
        } else {
            $consensus[$strand][$read][SEQ] = getSWConsensus(
                \@sortedSeq,
                $seqs,
                scalar(@{$downsample[$strand]})
            );
        }
        $consensus[$strand][$read][QUAL] = getConsensusQual_38($consensus[$strand][$read][SEQ]); 
    }    
}

# make a duplex consensus from two strand consensuses
sub makeDuplexConsensus {
    my (@reads) = @_;
    foreach my $read(@reads){
        if($consensus[STRAND1][$read][SEQ] and
           $consensus[STRAND2][$read][SEQ]){ # could be false if SSCS failed on a strand
            if($consensus[STRAND1][$read][SEQ] eq $consensus[STRAND2][$read][SEQ]){
                $consensus[DUPLEX][$read][SEQ] = $consensus[STRAND1][$read][SEQ];  
            } else {
                $consensus[DUPLEX][$read][SEQ] = maskTwoSequences($consensus[STRAND1][$read][SEQ], $consensus[STRAND2][$read][SEQ]);
            }
            $consensus[DUPLEX][$read][QUAL] = getConsensusQual_40($consensus[DUPLEX][$read][SEQ]);      
        } 
    }
}

# introduce IUPAC codes in place of base ambiguities to create a consensus between exactly two sequences
# the two sequences might be two input raw reads, or two previosly constructed merges or consensuses
# no identical sequences ever arrive at this sub, they were handled previously
sub maskTwoSequences {
    my ($ref, $qry) = @_; 

    # run fast Smith-Waterman on the sequence pair
    my ($qryOnRef, $score, $startQry, $endQry, $startRef) = smith_waterman($qry, $ref, 1);
    $score >= length($ref) * $minScoreFactor or return; # sequences did not match as well as required

    # M operations replace base conflicts with IUPAC codes for the pair of base values
    # D or I operations relative to reference are masked to N    
    # overruns of one sequence past another would be N bases and are not reported, i.e. are trimmed away
    foreach my $i(0..$#$qryOnRef){ 
        my $refBase = substr($ref, $i + $startRef, 1);
        $$qryOnRef[$i] ne $refBase and $$qryOnRef[$i] = $mOperations{$$qryOnRef[$i].$refBase} || 'N'; 
    }
    return join("", @$qryOnRef);
}

# use Smith-Waterman alignment against the most frequent sequence to resolve a consensus
# inputs to this sub are always a set of 3 or more input sequence values (possibly merged and thus with IUPAC codes)
sub getSWConsensus {
    my ($seqs, $quals, $nReadPairs) = @_;
    my $refLen = length($$seqs[0]);
    my $minScore = $refLen * $minScoreFactor;
    my $endRefExp = $refLen - 1;

    # run Smith-Waterman on all lesser-count sequences against the highest-count sequence
    my @sw = ([split("", $$seqs[0])]); # initialize with the reference sequence bases
    foreach my $i(1..$#$seqs){
        my ($qryOnRef, $score, $startQry, $endQry, $startRef, $endRef) = smith_waterman($$seqs[$i], $$seqs[0], 1);    
        $score >= $minScore or return; # fail SW on ANY read with too low a score
        $startRef > 0 and unshift @$qryOnRef, ("-") x $startRef; # ensure that all qry arrays are same length as ref (and thus as each other)
        $endRef < $endRefExp and push @$qryOnRef, ("-") x ($endRefExp - $endRef);
        push @sw, $qryOnRef;
    }
    
    # replace base conflicts with majority base value, or N if no consensus
    my $consensus;
    foreach my $posI(0..$endRefExp){ # cycle every position on the reference sequence
        my %baseCounts;
        foreach my $seqI(0..$#sw){   # cycle every sequence, including the reference
            $baseCounts{$sw[$seqI][$posI]} += @{$$quals{$$seqs[$seqI]}}; # account for all reads that had this same sequence
        }
        my @sortedBases = sort {$baseCounts{$b} <=> $baseCounts{$a}} keys %baseCounts;
        $consensus .= $baseCounts{$sortedBases[0]} >= $nReadPairs * $consensusFactor ?
            $sortedBases[0] :
            ($mOperations{join("", @sortedBases)} || 'N');  
    }
    $consensus =~ s/-//g; # strip away consensus D operation placeholders
    return $consensus;
}

1;

