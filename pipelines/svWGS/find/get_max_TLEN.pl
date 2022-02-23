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
    QUAL => 10,
    #-------------
    _PAIRED => 1,
    _PROPER => 2 # SAM FLAG bits
};

# working variables
my ($maxTlen, $prevQname, @alns) = (0);

# run the requested number of molecules
while (my $line = <STDIN>){
    my @f = split("\t", $line, SEQ + 2);
    if($prevQname and $f[QNAME] ne $prevQname){
        my $tLen = 0;
        if($alns[0][FLAG] & _PAIRED and
           $alns[0][FLAG] & _PROPER and 
           @alns == 2){ # unmerged, no supplementary, aligner has valid TLEN info
            $tLen = abs($alns[0][TLEN]);
        } elsif(@alns == 1) { # merged, TLEN is molecule length if not split
            $tLen = length($alns[0][SEQ]);
        }
        $maxTlen >= $tLen or $maxTlen = $tLen;
        @alns = ();
    }
    push @alns, \@f;
    $prevQname = $f[QNAME]; 
}

# return the result
print $maxTlen;
