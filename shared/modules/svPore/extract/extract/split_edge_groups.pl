use strict;
use warnings;

# separate training and SV edges into separate files

# load dependencies
our $script = "split_edge_groups";
our $error  = "$script error";
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
resetCountFile();

# environment variables
fillEnvVar(\our $N_CPU,  'N_CPU');
fillEnvVar(\our $EDGES_TMP_FILE, 'EDGES_TMP_FILE'); # the non-SV adapter training set
fillEnvVar(\our $EDGES_SV_FILE,  'EDGES_SV_FILE');

# open output handles
open my $tmpH,  "|-", "pigz -p $N_CPU -c | slurp -s 10M -o $EDGES_TMP_FILE" or die "could not open: $EDGES_TMP_FILE\n";
open my $svH,   "|-", "pigz -p $N_CPU -c | slurp -s 10M -o $EDGES_SV_FILE"  or die "could not open: $EDGES_SV_FILE\n";

# process data
my ($nReads, $nSv, $nNoSv, $prevQName, @edges) = (0, 0, 0);
while(my $edge = <STDIN>){
    chomp $edge;
    my ($qName) = split("\t", $edge, 2);     
    if($prevQName and $prevQName ne $qName){
        printMoleculeEdges();
        @edges = ();
    }
    push @edges, $edge;
    $prevQName = $qName;
}
printMoleculeEdges();
close $svH;
close $tmpH;

# print summary information
printCount($nReads, 'nReads',   'total reads processed');
printCount($nSv,    'nSv',      'reads with at least one candidate SV');
printCount($nNoSv,  'nNoSv',    'single-alignment reads with no SV');

# print a molecule's edges to the appropriate file(s)
sub printMoleculeEdges {
    $nReads++;
    if(@edges == 1){
        $nNoSv++;
        $nNoSv <= 10000 or return;
        print $tmpH join("\t", $edges[0], $nNoSv), "\n";
    } else {
        $nSv++;
        print $svH join("\n", map { join("\t", $_, $nSv) } @edges), "\n";
    }
}
