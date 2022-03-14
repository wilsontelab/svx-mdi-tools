use strict;
use warnings;

# initialize reporting
our $script = "coverage";

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(targets);
resetCountFile();

# environment variables
fillEnvVar(\my $MIN_MAPQ,        'MIN_MAPQ');
fillEnvVar(\my $GENOME_SIZE,     'GENOME_SIZE');
fillEnvVar(\our $TARGETS_BED,    'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING, 'REGION_PADDING');
fillEnvVar(\our $TARGET_SCALAR,  'TARGET_SCALAR', 1, 10); # use 10 bp target resolution for svCapture targets

# constants
use constant {
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    QUAL => 10,
    #-------------
    _PROPER => 2, # SAM FLAG bits
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
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
loadTargetRegions();
use vars qw($sumTargetLens $sumPaddedTargetLens %targetRegions);

# operating parameters
my $offTargetSize = $GENOME_SIZE - $sumPaddedTargetLens; # padded bases are not considered in coverage

# working variables
my ($prevName, $maxMapQ, @alns, @mol, %sumCoverage) = ("", 0);

# get TLEN from proper molecules
# don't multi-thread since we seek to sum over all source molecules
while(my $line = <STDIN>){
    
    # repeat line to STDOUT for saving BAM to disk
    print $line;
    
    # ignore SAM header lines
    $line =~ m/^@/ and next;

    # parse output one source molecule at a time
    my @aln = (split("\t", $line, TLEN + 2))[QNAME..TLEN];    
    if($prevName and $aln[QNAME] ne $prevName){
 
        # proceed if ANY alignment segment was of high enough MAPQ
        if($maxMapQ >= $MIN_MAPQ){

            # parse the information tags encoded in the read name
            @mol = split(":", $prevName);

            # identify pairs as proper and act accordingly
            if($mol[IS_MERGED]){ # merged reads, expect just one alignment; any supplemental = SV
                if(@alns == 1 and !($alns[0][FLAG] & _UNMAPPED)){
                    addProperMoleculeCoverage(getEnd($alns[0][POS], $alns[0][CIGAR]));
                } 
            } else { # unmerged reads, expect 2 alignments flagged as proper
                if(@alns == 2 and ($alns[0][FLAG] & _PROPER)){
                    addProperMoleculeCoverage(getEnd($alns[1][POS], $alns[1][CIGAR]));
                } 
            }         
        }

        # prepare for next read-pair
        ($maxMapQ, @alns) = (0);
    }

    # add new alignment to growing source molecule
    push @alns, \@aln;
    $maxMapQ >= $aln[MAPQ] or $maxMapQ = $aln[MAPQ];
    $prevName = $aln[QNAME];    
}

# print summary information
my $onTarget   = $sumTargetLens ? int($sumCoverage{'TT'} / $sumTargetLens * 1000 + 0.5) / 1000 : 0;
my $offTarget  = $offTargetSize ? int($sumCoverage{'--'} / $offTargetSize * 1000 + 0.5) / 1000 : 0;
my $enrichment = $offTarget     ? int($onTarget / $offTarget* 10 + 0.5) / 10 : 0;
printCount($sumCoverage{'TT'}, 'sumCoverage[TT]',  'number of bases covered in unpadded target regions');
printCount($sumCoverage{'--'}, 'sumCoverage[--]',  'number of bases covered in off-target regions');
printCount($onTarget,          'foldCoverage[TT]', 'coverage depth in target regions');
printCount($offTarget,         'foldCoverage[--]', 'coverage depth in off-target regions');
printCount($enrichment,        'enrichment',       'foldCoverage[TT] / foldCoverage[--]');

# add a proper molecule to the appropriate coverage sum
sub addProperMoleculeCoverage {
    my ($end) = @_;
    my $type = getTargetType($alns[0][POS], $end);
    $type and $sumCoverage{$type} += $end - $alns[0][POS] + 1;
    sub getTargetType {
        my ($pos1, $pos2) = @_;
        my $ct  = $targetRegions{$alns[0][RNAME]} or return "--";
        my $ct1 = $$ct{int($pos1 / $TARGET_SCALAR)};
        my $ct2 = $$ct{int($pos2 / $TARGET_SCALAR)};
        my $type1 = $ct1 ? $$ct1[1] : '';
        my $type2 = $ct2 ? $$ct2[1] : '';        
        if($type1 eq ON_TARGET or $type2 eq ON_TARGET){
            return 'TT'; # will include TA proper molecules crossing region edges
        } elsif($type1 or $type2){
            return "AA"; # proper molecules in padded regions will be ignored
        } else {
            return "--";
        }
    }    
}
