use strict;
use warnings;

# separate training and SV edges into separate files
# interleave SV edge groups for each read with a fake separator junction to aid downstream calculations in R

# load dependencies
our $script = "split_edge_groups";
our $error  = "$script error";
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
resetCountFile();

# constants
use constant {
    QNAME => 0, # edge format fields
    NODE1 => 1,
    QSTART => 2,    
    NODE2 => 3,
    QEND => 4,
    MAPQ => 5,
    CIGAR => 6,
    GAP_COMPRESSED_IDENTITY => 7,
    EDGE_TYPE => 8,
    EVENT_SIZE => 9,
    INSERT_SIZE => 10,
    N_STRANDS => 11,
    #-------------
    baseQual => 12,   
    alnBaseQual => 13, 
    alnSize => 14, # added to edges by processRead_
    sStart => 15,
    sEnd => 16,
    # #-------------
    # clip5 => 0, # adapter scores added to edges here by addAdaptersScores
    # score5 => 0,
    # nBases5 => 0,
    # start5 => 0,
    # end5 => 0,
    # clip3 => 0,
    # score3 => 0,
    # nBases3 => 0,
    # start3 => 0,
    # end3 => 0,
    # score5C => 0,
    # nBases5C => 0,
    # start5C => 0,
    # end5C => 0,
    # score3C => 0,
    # nBases3C => 0,
    # start3C => 0,
    # end3C => 0,
    # #-------------
    # channel => 0, # added to edges by finishEdge
    # pod5File => 0,
    # blockN => 0,
    # edgeN => 0,
};

# environment variables
fillEnvVar(\our $EDGES_TMP_FILE, 'EDGES_TMP_FILE'); # the non-SV adapter training set
fillEnvVar(\our $EDGES_SV_FILE,  'EDGES_SV_FILE');

# open output handles
open my $tmpH,  "|-", "gzip -c | slurp -s 10M -o $EDGES_TMP_FILE" or die "could not open: $EDGES_TMP_FILE\n";
open my $svH,   "|-", "gzip -c | slurp -s 10M -o $EDGES_SV_FILE"  or die "could not open: $EDGES_SV_FILE\n";

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
