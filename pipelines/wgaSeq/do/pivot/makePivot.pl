use strict;
use warnings;

# constants
use constant {
	CHROM => 0, # BED fields
	START => 1,
	END_ => 2,
	BIN_N => 3, # is "." in input
	GC => 4,
	STRAND => 5,
	EXCLUDED => 6,
	GAP => 7,
	BAD_REGION => 8,
	UMAP => 9,
	GENMAP => 10
};

# arguments
my ($chrom) = @ARGV;

# working variables
my ($prevBin, %counts);
my @binColumns = qw(chrom start end binN gc strand excluded gap badRegion umap genmap);
my @cells = split(" ", $ENV{CELL_IDS});

# initialize output
$chrom eq 'chr1' and print join("\t", @binColumns, @cells), "\n";

# create a crosstab of genome bin x cell, values are unique read pair counts
while(my $line = <STDIN>){
	chomp $line;
	my @f = split("\t", $line);
	if ($prevBin and ($$prevBin[CHROM] ne $f[CHROM] or
					  $$prevBin[START] != $f[START])) {
		commitBin();
		%counts = ();
    }
    $counts{$f[$#f-3]}++; # unique read counts by cell in each bin
	$prevBin = \@f;
}
commitBin();

sub commitBin {
	print join("\t", @$prevBin[CHROM..GENMAP]);
	foreach my $cell(@cells){
		my $count = $counts{$cell} || 0;
		print "\t$count";
	}
	print "\n";
}
