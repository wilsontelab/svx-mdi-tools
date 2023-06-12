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
    # _QSTART => 2,    
    NODE2 => 3,
    # _QEND => 4,
    # _MAPQ => 5,
    # CIGAR => 6,
    # GAP_COMPRESSED_IDENTITY => 7,
    EDGE_TYPE => 8,
    # EVENT_SIZE => 9,
    # INSERT_SIZE => 10,
    # N_STRANDS => 11,
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
fillEnvVar(\our $PIPELINE_DIR,     'PIPELINE_DIR');
fillEnvVar(\our $WINDOW_SIZE,      'WINDOW_SIZE');
fillEnvVar(\our $GENOME_FASTA,     'GENOME_FASTA');
fillEnvVar(\our $EDGES_NO_SV_FILE, 'EDGES_NO_SV_FILE');
fillEnvVar(\our $QNAMES_FILE,      'QNAMES_FILE');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svPore";
map { require "$PIPELINE_DIR/extract/$_.pl" } qw(initialize_windows);
initializeWindowCoverage();

# open output handles
open my $nosvH,    "|-", "gzip -c | slurp -s 10M -o $EDGES_NO_SV_FILE" or die "could not open: $EDGES_NO_SV_FILE\n";
open my $qNamesH,  "|-", "slurp -s 10M -o $QNAMES_FILE" or die "could not open: $QNAMES_FILE\n";

# process data
my ($nTmpLines, $prevQName, @lines) = (0);
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line, 11);
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
close $qNamesH;

# print a molecule's edges and/or qNames to the appropriate file(s)
sub printMoleculeEdges {
    if(@lines == 1){
        my $line = join("\t", @{$lines[0]})."\n";
        print $nosvH $line; # the record of all simple alignment edges at low base resolution, without CIGAR strings
        if($nTmpLines < 10000){
            print $qNamesH $prevQName."\n"; # qNames of 10K non-SV reads for training adapter models
            $nTmpLines++;
        }
    } else {
        print $qNamesH $prevQName."\n"; # qNames of all initial candidate SV reads
    }
}
