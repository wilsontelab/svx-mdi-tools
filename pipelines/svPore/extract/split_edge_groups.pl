use strict;
use warnings;

# separate training and SV edges into separate files
# interleave SV edge groups for each read with a fake separator junction to aid downstream calculations in R

my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);

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
    passedFlanks => 12, # added to edges by processRead_
    baseQual => 13,
    sStart => 14,
    sEnd => 15,
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
my ($prevQName, @edges);
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

# print a molecule's edges to the appropriate file(s)
my $edgeSetI = 0;
# my @spacer;
sub printMoleculeEdges {
    if(@edges == 1){
        print $tmpH join("\t", $edges[0], 0), "\n";
    } else {
        # if($edgeSetI){
        #     $edgeSetI++;
        #     print $svH join("\t", @spacer, $edgeSetI), "\n";
        # } else {
        #     my @x = split("\t", $edges[0]);
        #     @spacer = ("NA") x scalar(@x);
        #     $spacer[EDGE_TYPE] = "S";
        #     $spacer[$#spacer - 1] = 1; # blockN
        #     $spacer[$#spacer] = 1;     # edgeN
        # }
        $edgeSetI++;
        print $svH join("\n", map { join("\t", $_, $edgeSetI) } @edges), "\n";
    }
}
