use strict;
use warnings;

# apply fuzzy positions to keep read pairs matching an amplicon index

use constant {
    CHROM1  => 0,
    SIDE1   => 1,
    POS1    => 2,
    CHROM2  => 3,
    SIDE2   => 4,
    POS2    => 5,
    #--------------
    COUNT   => 6,
    AMPLICON=> 7,
    INDEX   => 8,
    #--------------
    UMI1    => 6,
    UMI2    => 7,
    SEQ1    => 8,
    QUAL1   => 9,
    SEQ2    => 10,
    QUAL2   => 11,
    OVERLAP => 12,
    MERGED  => 13
};

# collect the list of allowed outer node pairs
my %allowedNodePairs;
open my $inH, "<", $ENV{ALLOWED_FILE} or die "not found: $ENV{ALLOWED_FILE}\n";
while(my $line = <$inH>){
    chomp $line;
    my @line = split("\t", $line);
    my $key = join("\t", @line[CHROM1..POS2]);
    $allowedNodePairs{$key} = $line[AMPLICON];
}
close $inH;

# filter the stream
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line);
    my $key = join("\t", @line[CHROM1..POS2]);
    $allowedNodePairs{$key} and print join("\t", 
        $allowedNodePairs{$key},     
        @line[SEQ1, SEQ2, QUAL1, QUAL2, OVERLAP]
    )."\n";
}
