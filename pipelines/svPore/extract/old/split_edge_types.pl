use strict;
use warnings;

# finalize the parsing of re-extracted edges (all SV, some no-sv for adapter training)
# by splitting them into separate files

# constants
use constant {
    QNAME => 0,
    # NODE1 => 1,
    # _QSTART => 2,    
    # NODE2 => 3,
    # _QEND => 4,
    # _MAPQ => 5,
    # CIGAR => 6,
    # GAP_COMPRESSED_IDENTITY => 7,
    # EDGE_TYPE => 8,
    # EVENT_SIZE => 9,
    # INSERT_SIZE => 10,
    # N_STRANDS => 11,
};

# open output handles
open my $nosvH,  "|-", "gzip -c | slurp -s 10M -o $ENV{EDGES_TMP_FILE}" or die "could not open: $ENV{EDGES_TMP_FILE}\n";
open my $svH,    "|-", "gzip -c | slurp -s 10M -o $ENV{EDGES_SV_FILE}"  or die "could not open: $ENV{EDGES_SV_FILE}\n";

# process data
my ($prevQName, @lines);
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line, 2);
    if($prevQName and $prevQName ne $line[QNAME]){
        printMoleculeEdges();
        @lines = ();
    }
    push @lines, $line;
    $prevQName = $line[QNAME];
}
printMoleculeEdges();
close $nosvH;
close $svH;

# print a molecule's edges to the appropriate file(s)
sub printMoleculeEdges {
    if(@lines == 1){
        print $nosvH $lines[0]."\n";
    } else {
        print $svH join("\n", @lines), "\n";
    }
}
