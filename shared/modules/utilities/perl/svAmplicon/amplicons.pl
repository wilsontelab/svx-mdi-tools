use strict;
use warnings;

# load 
use vars qw($COLLATE_PREFIX);
our @amplicons;
open my $inH, "<", "$COLLATE_PREFIX.amplicons.txt" or die "missing file: $COLLATE_PREFIX.amplicons.txt\n";
while(my $line = <$inH>){
    chomp $line;
    my @line = split("\t", $line);
    $amplicons[$line[0]] = \@line;
}
close $inH;
