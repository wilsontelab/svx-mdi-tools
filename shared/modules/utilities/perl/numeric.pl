use strict;
use warnings;

# utility functions that manipulate numbers

#----------------------------------------------------------
# number formatting
#----------------------------------------------------------
sub commify {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}
sub roundCount {
    my ($count, $scalar) = @_;
    $scalar or $scalar = 1000;
    int($count * $scalar + 0.5) / $scalar;
}
sub roundCount2 {
    my ($count) = @_;
    int($count * 100 + 0.5) / 100;
}

#----------------------------------------------------------
# aggregate functions
#----------------------------------------------------------
sub min {
    my ($v1, $v2) = @_;
    $v1 <= $v2 ? $v1 : $v2;
}
sub max {
    my ($v1, $v2) = @_;
    $v1 >= $v2 ? $v1 : $v2;
}
sub median {
    my (@data) = sort {$a <=> $b} @_;
    my $i = @data / 2;
    my $upper = $data[$i];
    @data % 2 and return $upper;
    my $lower = $data[$i - 1];
    return($lower + ($upper - $lower) / 2);
}
sub mean{
    my (@data) = @_;
    @data == 0 and die "no values sent to mean\n";
    my $sum = 0;
    foreach (@data) { $sum += $_ }
    return $sum / @data;
}
sub stdev{
    my (@data) = @_;
    @data == 1 and return 0;
    my $mean = mean(@data);
    my $sqSum = 0;
    foreach(@data) { $sqSum += ($mean-$_) ** 2 }
    return ($mean, ($sqSum / (@data-1)) ** 0.5);
}

1;
