use strict;
use warnings;

# action:
#     create a map of base-level coverage in a genome
# input format:
#     sorted BED3 plus MOL_CLASS
# output:
#     baseCoverage.breaks   = 0-referenced positions where coverage changes
#     baseCoverage.counts   = the corresponding coverage at each break
#     baseCoverage.index.gz = index file for extracting the coverage map of a chunk of a chromosome

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
require "$perlUtilDir/numeric.pl";

# constants
use constant {
	CHROM_INDEX => 0,
	START_POS_0 => 1, # i.e, 0-indexed start pos
	END_POS_1   => 2, # i.e, 1-indexed end pos
	MOL_CLASS   => 3, # either P (proper) or V (variant)
	#---------------------
	BUFFER_LEN   => 10000, # number of base positions in the circular buffer (10K is longer than ever needed for short-read pairs)
	CHUNK_SIZE   => 65536, # size of a chromosome span with progressive numbering == 65536 bp, to control map file size
	MAX_COVERAGE => 255,   # largest allowed fragment coverage, to control map file size
};

# operating variables
my @breaks = (0) x BUFFER_LEN; # positions where the genome base coverage changes, RESETS EVERY CHROM
my $nBreaksPrinted = 0; # number of breaks printed to coverage map, summed over all chroms
my $chunkIndex = 0;     # index of the working chrom chunk, RESETS EVERY CHROM
my $maxPosInChunk = CHUNK_SIZE - 1; # the last position in the working chrom chunk, RESETS EVERY CHROM
my $coverage = 0; # running coverage value as breaks are handled, RESETS EVERY CHROM
my $pos0 = -1;    # 0-indexed base position for first element of @breaks, RESETS EVERY CHROM
my $maxEnd1 = -1; # rightmost span end encountered so far, RESETS EVERY CHROM
my $nProperSpans = 0; # for expressing insertSizes as frequencies, summed over all chroms
my ($prevChromI, @insertSizes); # tally of TLENs for histogram

# output files
my $filePrefix = "$ENV{DATA_FILE_PREFIX}.baseCoverage";
my $posFile = "$filePrefix.breaks";
my $covFile = "$filePrefix.counts";
my $idxFile = "$filePrefix.index.gz";
open my $posH, ">", $posFile or die "could not open $posFile for writing\n$!\n";
open my $covH, ">", $covFile or die "could not open $covFile for writing\n$!\n";
open my $idxH, "|-", "gzip -c > $idxFile" or die "could not open $idxFile for writing\n$!\n";
print $idxH join("\t", qw(chromIndex chunkIndex cumNBreaks)), "\n";

# run all sorted spans
while (<STDIN>) {
	chomp;
	my ($chromI, $start0, $end1, $molClass) = split("\t");

	# count proper molecule insert sizes
	if($molClass eq 'P'){
		$nProperSpans++;
		$insertSizes[$end1 - $start0]++;
	} 

	# finish a chromosome
	if($prevChromI){
		if($prevChromI != $chromI){
			commitBreaks($maxEnd1 + 1); # commit any/all remaining breaks on the chromosome
			print $idxH join("\t", $prevChromI, $chunkIndex, $nBreaksPrinted), "\n";
			@breaks = (0) x BUFFER_LEN;
			$chunkIndex = 0;
			$maxPosInChunk = CHUNK_SIZE - 1;
			$coverage = 0;
			$pos0 = -1;
			$maxEnd1 = -1;

		# commit any break positions prior to the first position of the current span
		} else {
			$start0 > $pos0 and commitBreaks($start0);
		}
	}

	# initialize $pos0 to the start of the first span on the chromosome
	$pos0 == -1 and $pos0 = $start0;

	# construct the coverage map from both proper molecules and SV alignment spans
	$breaks[$start0 % BUFFER_LEN] += 1; # 0-referenced start of this mapped span, increment coverage up
	$breaks[$end1   % BUFFER_LEN] -= 1; # 0-referenced base position _after_ this span ends, decrement coverage down

	# prepare for next span
	$prevChromI = $chromI;
	$maxEnd1 > $end1 or $maxEnd1 = $end1;
}
commitBreaks($maxEnd1 + 1); # finish the last chromosome
print $idxH join("\t", $prevChromI, $chunkIndex, $nBreaksPrinted), "\n";
close $posH;
close $covH;
close $idxH;
printInsertSizes();

# parse and print the completed coverage spans
sub commitBreaks {
    my ($start0) = @_;
	foreach my $p0($pos0..min($maxEnd1, $start0 - 1)){
		my $increment = $breaks[$p0 % BUFFER_LEN] or next;
		while($p0 > $maxPosInChunk){
			print $idxH join("\t", $prevChromI, $chunkIndex, $nBreaksPrinted), "\n";
			$chunkIndex++;
			$maxPosInChunk = ($chunkIndex + 1) * CHUNK_SIZE - 1;
		}	
		$coverage += $increment;
		print $posH pack("S", $p0); # implicitly applies modulo to pos
		print $covH pack("C", $coverage > MAX_COVERAGE ? MAX_COVERAGE : $coverage);
		$nBreaksPrinted++;
		$breaks[$p0 % BUFFER_LEN] = 0;
	}
	$pos0 = $start0;                          
}

# print the insert size histogram
sub printInsertSizes {
	my $outFile = "$ENV{EXTRACT_PREFIX}.insertSizes.txt";
	open my $outH, ">", "$outFile" or die "\ncould not open file for writing:\n$outFile:\n$!\n";
	print $outH join("\t", 'insertSize', 'frequency'), "\n";
	foreach my $insertSize(0..$#insertSizes){
		my $count = $insertSizes[$insertSize] || 0;
		print $outH join("\t", $insertSize, $count / $nProperSpans), "\n";
	}
	close $outH;
}

1;
