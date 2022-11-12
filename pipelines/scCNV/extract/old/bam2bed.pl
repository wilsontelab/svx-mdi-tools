use strict;
use warnings;

# prepare alignments in bam file for pivot by parsing to BED
# each molecule is counted in the bin corresponding to the start of the molecule span
# this creates a trivial uncertainty at the edges of very large bins

# constants
use constant {
    QNAME => 0, # sam fields
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
	TAGS => 11
};

# bam2bed
while(my $line = <STDIN>){
	my @f = split("\t", $line, 12);
	$f[TAGS] =~ m/CB:Z:(\S+)/;
    print join("\t", $f[RNAME], $f[POS] - 1, $f[POS], $1, 0, '.'), "\n";
}
