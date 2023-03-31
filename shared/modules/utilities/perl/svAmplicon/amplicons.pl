use strict;
use warnings;

# load 
use vars qw($EXTRACT_PREFIX);
our @amplicons;
open my $inH, "<", "$EXTRACT_PREFIX.amplicons.txt" or die "missing file: $EXTRACT_PREFIX.amplicons.txt\n";
while(my $line = <$inH>){
    chomp $line;
    my @line = split("\t", $line);
    $amplicons[$line[0]] = \@line;
}
close $inH;
