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
    QNAME => 0,
    NODE1 => 1,
    CIGAR1 => 2,
    BLAST_IDENTITY1 => 3,
    GAP_COMPRESSED_IDENTITY1 => 4,
    NODE2 => 5,
    CIGAR2 => 6,
    BLAST_IDENTITY2 => 7,
    GAP_COMPRESSED_IDENTITY2 => 8,
    EDGE_TYPE => 9,
    # _MAPQ => 10,
    # SV_SIZE => 11,
    # INSERT_SIZE => 12,
    # XSTART => 13,
    # XEND => 14,
    # EDGE_CLASS => 15,
    # N_STRANDS => 16,
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
fillEnvVar(\our $EDGES_NO_SV_FILE, 'EDGES_NO_SV_FILE');
fillEnvVar(\our $EDGES_SV_FILE,    'EDGES_SV_FILE');
fillEnvVar(\our $EDGES_TMP_FILE,   'EDGES_TMP_FILE');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svPore";
map { require "$ACTION_DIR/extract/$_.pl" } qw(initialize_windows);
initializeWindowCoverage();

# open output handles
open my $nosvH, "|-", "gzip -c | slurp -s 10M -o $EDGES_NO_SV_FILE" or die "could not open: $EDGES_NO_SV_FILE\n";
open my $svH,   "|-", "gzip -c | slurp -s 10M -o $EDGES_SV_FILE"    or die "could not open: $EDGES_SV_FILE\n";
open my $tmpH,  "|-", "gzip -c | slurp -s 10M -o $EDGES_TMP_FILE"   or die "could not open: $EDGES_TMP_FILE\n";

# process data
my ($nTmpLines, $prevQName, @lines) = (0);
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line, 11);
    if($line[EDGE_TYPE] ne ALIGNMENT){
        $line[CIGAR1] = "NA";
        $line[CIGAR2] = "NA";
    }
    if($prevQName and $prevQName ne $line[QNAME]){
        printMoleculeEdges();
        @lines = ();
    }
    push @lines, \@line;
    $prevQName = $line[QNAME];
    $line[EDGE_TYPE] eq ALIGNMENT or next;
    incrementWindowCoverage(@line[NODE1, NODE2]);
}
printMoleculeEdges();
printWindowCoverage();
close $nosvH;
close $svH;
close $tmpH;

# print a molecule's edges to the appopriate file(s)
sub printMoleculeEdges {
    if(@lines == 1){
        my $line = join("\t", @{$lines[0]})."\n";
        print $nosvH $line;
        if($nTmpLines < 10000){
            print $tmpH $line;
            $nTmpLines++;
        }
    } else {
        print $svH join("\n", map { join("\t", @$_) } @lines), "\n";
    }
}
