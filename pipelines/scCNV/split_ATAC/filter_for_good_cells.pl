
# load the good cells
my %goodCells;
open my $inH, "<", $ENV{GOOD_CELLS_FILE} or die "could not open: $ENV{GOOD_CELLS_FILE}\n";
while (my $line = <$inH>){
    chomp $line;
    my @f = split(",", $line);
    $goodCells{$f[0]} = $f[1];
}
close $inH;

# filter alignments in the bam stream
while (my $line = <STDIN>){
    $line =~ m/CB:Z:(\S+)/ or next;
    $goodCells{$1} or next;
    print join("\t", $goodCells{$1}, $line);
}
