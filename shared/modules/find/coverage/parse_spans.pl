use strict;
use warnings;

# action:
#     expand de-duplicated molecules to spans for coverage map
# output format:
#     BED3 plus MOL_CLASS

use constant {
    MOL_CLASS => 0,
    MOL_STRAND => 1,
    NODE1 => 2,
    NODE2 => 3,
    SPANS => 4,
    #------------
    UMI => 0,
    RNAME_INDEX => 1,
    _SIDE => 2,
    _POS => 3,
    IS_CLIPPED => 4,
};

while(<STDIN>){
    chomp;
    my @F = split("\t");

    # for proper molecules use the entire span of each source molecule
    if($F[SPANS] eq '-'){
        my @node1 = split(":", $F[NODE1]);
        my @node2 = split(":", $F[NODE2]);
        print join("\t", $node1[RNAME_INDEX], $node1[_POS] - 1, $node2[_POS], 'P'), "\n";

    } elsif($F[MOL_CLASS] eq 'P') {
        my @node1 = split(":", $F[NODE1], 3);        
        my @pos = split(":", $F[SPANS]);
        print join("\t", $node1[RNAME_INDEX], $pos[0] - 1, $pos[1], 'P'), "\n";

    # for SV molecules, use the span of all alignments, but not the gaps between them
    } else {
        $F[SPANS] =~ s/::/\tV\n/g;
        $F[SPANS] =~ s/:/\t/g;
        print $F[SPANS], "\tV\n";
    }
}
