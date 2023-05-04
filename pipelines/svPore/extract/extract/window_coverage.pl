use strict;
use warnings;

# finalize the genome fragment coverage map

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# constants
use constant {
    QNAME => 0, # PAF fields
    NODE1 => 1,
    CIGAR1 => 2,
    NODE2 => 3,
    CIGAR2 => 4,
    EDGE_TYPE => 5,
    # _MAPQ => 6,
    # SV_SIZE => 7,
    # INSERT_SIZE => 8,
    # XSTART => 9,
    # XEND => 10,
    # EDGE_CLASS => 11,
    # N_STRANDS => 12,
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I"
};

# environment variables
fillEnvVar(\our $EXTRACT_PREFIX,   'EXTRACT_PREFIX');
fillEnvVar(\our $ACTION_DIR,       'ACTION_DIR');
fillEnvVar(\our $N_CPU,            'N_CPU'); # user options, or derived from them
fillEnvVar(\our $WINDOW_SIZE,      'WINDOW_SIZE');
fillEnvVar(\our $MIN_SV_SIZE,      'MIN_SV_SIZE');
fillEnvVar(\our $GENOME_FASTA,     'GENOME_FASTA');
fillEnvVar(\our $USE_CHR_M,        'USE_CHR_M');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svPore";
map { require "$ACTION_DIR/extract/$_.pl" } qw(initialize_windows);
initializeWindowCoverage();

# process data
while(my $line = <STDIN>){
    print $line; # repeat for nodes file
    chomp $line;
    my @line = split("\t", $line);
    $line[EDGE_TYPE] eq ALIGNMENT or next;
    incrementWindowCoverage(@line[NODE1, NODE2]);
}
printWindowCoverage();
