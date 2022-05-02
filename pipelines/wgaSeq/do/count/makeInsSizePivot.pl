use strict;
use warnings;

# working variables
my @cells = split(" ", $ENV{CELL_IDS});
my @counts;

# initialize output
print join("\t", 'insSize', @cells), "\n";

# create a crosstab of insSize x cell, values are unique read pair counts
while(my $line = <STDIN>){
	chomp $line;
	my ($insSize, $cell) = split("\t", $line);
	$insSize > 0 or next; # failsafe
	$counts[$insSize]{$cell}++;
	$counts[0]{$cell}++;
}

# print the final table
foreach my $insSize(1..1500){
	my $c = $counts[$insSize];
	$c or $c = {};
	print $insSize;
	foreach my $cell(@cells){
		my $count = $$c{$cell} || 0;
		my $N = $counts[0]{$cell} || 0;
		my $freq = $N ? $count / $N : 0;
		print "\t$freq";
	}
	print "\n";
}
