use strict;
use warnings;

use constant {
    AMP_AMPLICON_ID => 0, # amplicon fields
    AMP_PROPER => 1,
    AMP_MOL_COUNT => 2,
    AMP_CHROM1 => 3,
    AMP_SIDE1 => 4,
    AMP_POS1 => 5,
    AMP_REF1 => 6,
    AMP_CHROM2 => 7,
    AMP_SIDE2 => 8,
    AMP_POS2 => 9,
    AMP_REF2 => 10,
};

# load 
use vars qw($COLLATE_PREFIX);
our @amplicons;
open my $inH, "<", "$COLLATE_PREFIX.amplicons.txt" or die "missing file: $COLLATE_PREFIX.amplicons.txt\n";
while(my $line = <$inH>){
    chomp $line;
    my @line = split("\t", $line);
    $amplicons[$line[AMP_AMPLICON_ID]] = \@line;
}
close $inH;
