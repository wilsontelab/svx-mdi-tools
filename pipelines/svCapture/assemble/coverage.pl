use strict;
use warnings;

# profile the coverage of padded target regions by proper molecules

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(targets);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# environment variables and arguments passed by assemble.R
our($GENOME, $TARGETS_BED, $REGION_PADDING, $MIN_MAPQ, $NAME_BAM_FILE) = @ARGV;
fillEnvVar(\our $MIN_FLANK_SIZE,    'MIN_FLANK_SIZE');
fillEnvVar(\our $MIN_FAMILY_SIZE,   'MIN_FAMILY_SIZE');
fillEnvVar(\our $TARGET_SCALAR,     'TARGET_SCALAR', 1, 10); # use 10 bp target resolution for matching svCapture targets ...
fillEnvVar(\our $COVERAGE_BIN_SIZE, 'COVERAGE_BIN_SIZE');    # ... but create coverage profiles at the requested bin size

# constants
use constant {
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    # RNEXT => 6,
    # PNEXT => 7,
    # TLEN => 8,
    # SEQ => 9,
    # QUAL => 10,
    #-------------
    _PROPER => 2, # SAM FLAG bits
    _UNMAPPED => 4,
    #-------------
    MOL_ID => 0, # @mol fields, carried here in QNAME
    UMI1 => 1,
    UMI2 => 2,
    IS_MERGED => 3,
    IS_DUPLEX => 4,
    STRAND_COUNT1 => 5,
    STRAND_COUNT2 => 6,
    #-------------
    ON_TARGET   => 'T',
    NEAR_TARGET => 'A', # A for adjacent
    OFF_TARGET  => '-',
};

# load the target regions
my $regions = getTargetRegions();
loadTargetRegions(1);
use vars qw($sumTargetLens %targetRegions);

# working variables
my ($prevName, $maxMapQ, @alns, @mol, %sums, @insertSizes, @coverage) = ("", 0);

# read bam stream
open my $inH, "-|", "samtools view $NAME_BAM_FILE | cut -f 1-6" or die "could not open alignment stream: $!\n";
while(my $line = <$inH>){

    # parse output one source molecule at a time
    chomp $line;
    my @aln = split("\t", $line);    
    if($prevName and $aln[QNAME] ne $prevName){
        processSourceMolecule();
        ($maxMapQ, @alns) = (0);
    }

    # add new alignment to growing source molecule
    push @alns, \@aln;
    $maxMapQ >= $aln[MAPQ] or $maxMapQ = $aln[MAPQ];
    $prevName = $aln[QNAME];    
}
processSourceMolecule();
close $inH;

# return sample aggregates
my @countColumns    = qw(n_source_molecules n_on_target n_on_target_filtered);
my @coverageColumns = qw(kbp_on_target kbp_on_target_effective kbp_on_target_filtered kbp_on_target_filtered_effective);
my @binColumns      = qw(bin_counts);
my $header = join("\t", 
    "sumTargetLens", 
    "N50",
    @countColumns,
    @coverageColumns,
    @binColumns
);
my $data = join("\t", 
    $sumTargetLens, 
    getN50(@insertSizes), 
    (map { $sums{$_} || 0 } @countColumns), 
    (map { 
        my $bp = $sums{$_} || 0;
        int($bp / 1000 + 0.5); # return kbp to prevent int32 overruns
    } @coverageColumns),
    join("::", map {
        my $maxBinI = int(($$regions[$_]{paddedEnd} - $$regions[$_]{paddedStart}) / $COVERAGE_BIN_SIZE);
        my $regId = $_ + 1;
        join(",", map {
            $coverage[$regId][$_] || 0;
        } 0..$maxBinI)
    } 0..$#{$regions})
);
print "$header\n$data\n";

# process a source molecule if ANY alignment segment was of high enough MAPQ
# only continue with proper molecules for coverage profiling
sub processSourceMolecule {
    $sums{n_source_molecules}++;
    $maxMapQ >= $MIN_MAPQ or return;  
    @mol = split(":", $prevName);
    if($mol[IS_MERGED]){ # merged reads, expect just one alignment; any supplemental = SV
        if(@alns == 1 and !($alns[0][FLAG] & _UNMAPPED)){
            addProperMolecule(getEnd($alns[0][POS], $alns[0][CIGAR]));
        } 
    } else { # unmerged reads, expect 2 alignments flagged as proper
        if(@alns == 2 and ($alns[0][FLAG] & _PROPER)){
            $alns[0][POS] > $alns[1][POS] and @alns = reverse @alns;
            addProperMolecule(getEnd($alns[1][POS], $alns[1][CIGAR]));
        } 
    }         
}

# add information about a proper molecule within a target region
sub isOnTarget {
    my ($pos1, $pos2) = @_;
    my $ct  = $targetRegions{$alns[0][RNAME]} or return {};
    my $ct1 = $$ct{int($pos1 / $TARGET_SCALAR)};
    my $ct2 = $$ct{int($pos2 / $TARGET_SCALAR)};
    my $regionId = $ct1 ? $$ct1[0] : undef; # 1-referenced index of row of targets bed
    $regionId or return {};
    my $type1    = $ct1 ? $$ct1[1] : '';
    my $type2    = $ct2 ? $$ct2[1] : '';
    {
        regionId => $regionId,
        padded   => ($type1 eq ON_TARGET or $type2 eq ON_TARGET or $type1 eq NEAR_TARGET or $type2 eq NEAR_TARGET),
        unpadded => ($type1 eq ON_TARGET or $type2 eq ON_TARGET) # will include TT and TA proper molecules crossing region edges but not AA or -- proper
    };
}
sub getTargetBinI {
    my ($regId, $pos) = @_;
    int(($pos - $$regions[$regId - 1]{paddedStart}) / $COVERAGE_BIN_SIZE);
}
sub addProperMolecule {
    my ($end) = @_;
    my $target = isOnTarget($alns[0][POS], $end);

    # add reads mapped to padded target regions to binned coverage profiles
    $$target{padded} or return;
    my $regId = $$target{regionId};
    for(my $binI = getTargetBinI($regId, $alns[0][POS]); $binI <= getTargetBinI($regId, $end); $binI++){
        $binI >= 0 and $coverage[$regId][$binI]++;
    }

    # count and sum reads within the unpadded capture targets as overall estimates of library depth for inter-sample comparison
    $$target{unpadded} or return;
    my $insertSize = $end - $alns[0][POS] + 1;
    my $effectiveSize = max(0, $insertSize - 2 * $MIN_FLANK_SIZE);
    push @insertSizes, $insertSize;

    $sums{n_on_target}++;
    $sums{kbp_on_target} += $insertSize;
    $sums{kbp_on_target_effective} += $effectiveSize; # a more accurate coverage estimate that accounts for need to flank SV junctions during calling

    my $familySize = $mol[STRAND_COUNT1] + $mol[STRAND_COUNT2];
    $familySize >= $MIN_FAMILY_SIZE or return;

    $sums{n_on_target_filtered}++; # here, don't include reads the couldn't pass the SV filter family size threshold
    $sums{kbp_on_target_filtered} += $insertSize;
    $sums{kbp_on_target_filtered_effective} += $effectiveSize;
}
