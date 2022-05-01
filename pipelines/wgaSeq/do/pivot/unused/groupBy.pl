use strict;
use warnings;

# get columns
my ($maxDataCol, $valCol) = @ARGV;
$maxDataCol--;
$valCol--;

# constants
use constant {
    CHROM => 0,
    START => 1,
    END_ => 2,
    BIN_N => 3,
};

# working variables
my ($prevBinN, $sum, $n, $x) = (0, 0, 0);

# average the score over a bin
while (my $line = <STDIN>) {
    chomp $line;
    my @f = split("\t", $line);
    if ($prevBinN and $f[BIN_N] != $prevBinN) {
        print join("\t", @$x[0..$maxDataCol], $sum / $n), "\n";
        ($sum, $n, $x) = (0, 0);
    }
    $prevBinN = $f[BIN_N];
    $sum += $f[$valCol];
    $n++;
    !$x and $x = \@f;
}
print join("\t", @$x[0..$maxDataCol], $sum / $n), "\n";

