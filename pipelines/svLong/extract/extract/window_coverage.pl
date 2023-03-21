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
    NODE2 => 2,
    EDGE_TYPE => 3,
    _MAPQ => 4,
    SV_SIZE => 5,
    INSERT_SIZE => 6,
    N_STRANDS => 7,
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
fillEnvVar(\our $WINDOW_POWER,     'WINDOW_POWER');
fillEnvVar(\our $WINDOW_SIZE,      'WINDOW_SIZE');
fillEnvVar(\our $MIN_SV_SIZE,      'MIN_SV_SIZE');
fillEnvVar(\our $GENOME_FASTA,     'GENOME_FASTA');
fillEnvVar(\our $USE_CHR_M,        'USE_CHR_M');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svLong";
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
