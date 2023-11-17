use strict;
use warnings;

# action:
#     create a map of base-level coverage of _short read pairs_ in a genome
# input format:
#     sorted BED3 plus MOL_CLASS
# output:
#     baseCoverage.breaks   = 0-referenced positions where coverage changes
#     baseCoverage.counts   = the corresponding coverage at each break start until the next break
#     baseCoverage.index.gz = index file for extracting the coverage map of a chunk of a chromosome
#                             also contains a low resolution bin/chunk-based coverage map

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
	BUFFER_LEN   => 10000,     # number of base positions in the circular buffer (10K is longer than ever needed for short-read pairs)
	CHUNK_SIZE   => 2**16,     # size of a chromosome span with progressive numbering == 65536 bp, to control map file size
	MAX_COVERAGE => 2**16 - 1, # largest allowed fragment coverage, to control map file size
};

# operating variables
my @breaks = (0) x BUFFER_LEN; # positions where the genome base coverage changes, RESETS EVERY CHROM
my $nBreaksPrinted = 0; # number of breaks printed to coverage map, summed over all chroms
my $chunkIndex = 0;     # index of the working chrom chunk, RESETS EVERY CHROM
my $maxPosInChunk = CHUNK_SIZE - 1; # the last position in the working chrom chunk, RESETS EVERY CHROM
my $coverage = 0; # running coverage value as breaks are handled, RESETS EVERY CHROM
my $pos0 = -1;    # 0-indexed base position for first element of @breaks, RESETS EVERY CHROM
my $maxEnd1 = -1; # rightmost span end encountered so far, RESETS EVERY CHROM
my $cumChunkWeightedCoverage = 0; # for calculating the mean coverage over all bases in each chunk
my $cumChunkWeights = 0;
my $prevBreakPos = 1;
my ($prevChromI); 

# output files
my $posFile = "$ENV{COVERAGE_PREFIX}.breaks";
my $covFile = "$ENV{COVERAGE_PREFIX}.counts";
my $idxFile = "$ENV{COVERAGE_PREFIX}.index.gz";
open my $posH, ">", $posFile or die "could not open $posFile for writing\n$!\n";
open my $covH, ">", $covFile or die "could not open $covFile for writing\n$!\n";
open my $idxH, "|-", "gzip -c > $idxFile" or die "could not open $idxFile for writing\n$!\n";
print $idxH join("\t", qw(chromIndex chunkIndex cumNBreaks coverage)), "\n";

# run all sorted spans
while (<STDIN>) {
	chomp;
	my ($chromI, $start0, $end1, $molClass) = split("\t"); 

	# finish a chromosome
	if($prevChromI){
		if($prevChromI != $chromI){
			commitBreaks($maxEnd1 + 1); # commit any/all remaining breaks on the chromosome
			print $idxH join("\t", $prevChromI, $chunkIndex, $nBreaksPrinted, getChunkCoverage()), "\n";
			@breaks = (0) x BUFFER_LEN;
			$chunkIndex = 0;
			$maxPosInChunk = CHUNK_SIZE - 1;
			$coverage = 0;
			$pos0 = -1;
			$maxEnd1 = -1;
			$prevBreakPos = 1;

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
print $idxH join("\t", $prevChromI, $chunkIndex, $nBreaksPrinted, getChunkCoverage()), "\n";
close $posH;
close $covH;
close $idxH;

# parse and print the completed coverage spans
sub commitBreaks {
    my ($start0) = @_;
	foreach my $p0($pos0..min($maxEnd1, $start0 - 1)){
		my $increment = $breaks[$p0 % BUFFER_LEN] or next;

		# handle chunk breaks
		while($p0 > $maxPosInChunk){
			my $breakWeight = $maxPosInChunk - $prevBreakPos + 1;
			$cumChunkWeightedCoverage += $coverage * $breakWeight;
			$cumChunkWeights += $breakWeight;
			$prevBreakPos = $maxPosInChunk + 1;	
			print $idxH join("\t", $prevChromI, $chunkIndex, $nBreaksPrinted, getChunkCoverage()), "\n";
			$chunkIndex++;
			$maxPosInChunk = ($chunkIndex + 1) * CHUNK_SIZE - 1;
		}

		# keep track of the mean coverage per chunk prior to this break
		my $breakWeight = $p0 - $prevBreakPos;
        $cumChunkWeightedCoverage += $coverage * $breakWeight;
		$cumChunkWeights += $breakWeight;
		$prevBreakPos = $p0;

		# increment and commit this break 
		$coverage += $increment;
		print $posH pack("S", $p0); # implicitly applies modulo to pos
		print $covH pack("S", $coverage > MAX_COVERAGE ? MAX_COVERAGE : $coverage);
		$nBreaksPrinted++;
		$breaks[$p0 % BUFFER_LEN] = 0;
	}
	$pos0 = $start0;                          
}
sub getChunkCoverage {
	my $chunkMeanCoverage = $cumChunkWeights ? $cumChunkWeightedCoverage / $cumChunkWeights : 0;
	$cumChunkWeightedCoverage = 0;
	$cumChunkWeights = 0;
	int($chunkMeanCoverage * 10 + 0.5) / 10; # record chunk coverage to 0.1 precision
}

1;
