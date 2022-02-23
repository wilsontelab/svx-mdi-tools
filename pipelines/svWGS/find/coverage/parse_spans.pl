use strict;
use warnings;

# action:
#     expand de-duplicated molecules to spans for coverage map and insert size analysis
# input format:
#     class   strand  chr1    side1   pos1            chr2    side2   pos2            spans
#     P       0       14      R       52363064        14      L       52363453        -
#     P       0       14      R       52363064        14      L       52363453        52363074:52363453
#     V       0       13      R       96113929        17      R       76065242        17:76065241:76065392::13:96113928:96114079
# output format:
#     BED3 plus MOL_CLASS

# constants
use constant {
    MOL_CLASS   => 0,
    STRAND  => 1,
    CHR1    => 2,
    SIDE1   => 3,
    POS1    => 4,
    CHR2    => 5,
    SIDE2   => 6,
    POS2    => 7,
    SPANS   => 8
};

while(<STDIN>){
    chomp;
    my @F = split("\t");

    # for proper molecules use the entire span of each source molecule
    if($F[SPANS] eq '-'){
        print join("\t", $F[CHR1], $F[POS1] - 1, $F[POS2], 'P'), "\n";

    } elsif($F[MOL_CLASS] eq 'P') {
        my @pos = split(":", $F[SPANS]);
        print join("\t", $F[CHR1], $pos[0] - 1, $pos[1], 'P'), "\n";

    # for SV molecules, use the span of all alignments, but not the gaps between them
    } else {
        $F[SPANS] =~ s/::/\tV\n/g;
        $F[SPANS] =~ s/:/\t/g;
        print $F[SPANS], "\tV\n";
    }
}
