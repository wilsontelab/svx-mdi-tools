use strict;
use warnings;

# STEP 1 - establish coherent merge levels of all molecules in the group

# we cannot easily build a consensus on a combination or merged + unmerged read pairs
# therefore, reject unmerged (presumed lower quality) if merged have same signature
# if only one strand was merged, demote from duplex to single-strand molecule

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

# working variables
use vars qw(@readPairs $mol @strands);
my $nullSymbol = '*';

sub parseMergeLevels {
    
    # initialize
    $$mol[IS_MERGED] = 0;
    
    # determine which strands have how many merged molecules
    my @merged = (0, 0);
    foreach my $strand(@strands){
        foreach my $readPair(@{$readPairs[$strand]}){
            if($$readPair[SEQ2] eq $nullSymbol){
                $merged[$strand]++;
                $$readPair[SEQ_MERGED]  = $$readPair[SEQ1];
                $$readPair[QUAL_MERGED] = $$readPair[QUAL1];
            } 
        }
        
        # if mix of merged and unmerged, discard the unmerged (presumed lower quality) read pairs
        if($merged[$strand]){
            $$mol[IS_MERGED] = 1;
            $merged[$strand] < @{$readPairs[$strand]} and # merged are always sorted first
                @{$readPairs[$strand]} = @{$readPairs[$strand]}[0..($merged[$strand]-1)];
        }
    }
    
    # if two strands but only one was merged, discard the unmerged strand
    if($$mol[IS_MERGED] and
       @strands == 2 and
       (!$merged[STRAND1] or !$merged[STRAND2])){
        @strands = $merged[STRAND1] ? (STRAND1) : (STRAND2);
    }
}

1;
