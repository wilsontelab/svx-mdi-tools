use strict;
use warnings;

# write a FASTQ stream of kept and grouped molecule sequences suitable for re-alignment
# each molecule may have one merged or two unmerged reads
# name = molId:ampliconId:nOverlapBases:molCount:merged:readN

use constant {
    AMPLICON=> 0,
    SEQ1    => 1,
    SEQ2    => 2,
    QUAL1   => 3,
    QUAL2   => 4,
    OVERLAP => 5,
    COUNT   => 6
};

# parse the stream
my $i = 0;
while(my $line = <STDIN>){
    $i++;
    chomp $line;
    my @line = split("\t", $line);
    my $name = join(":", $i, @line[AMPLICON, OVERLAP, COUNT], $line[SEQ2] eq "*" ? 1 : 0);
    my $qual1 = $line[COUNT] > 1 ? getModalQual($line[QUAL1]) : $line[QUAL1]; # repeat sequences assumed to have high base quality
    print "\@$name:1\n$line[SEQ1]\n+\n$qual1\n";
    if($line[SEQ2] ne "*"){
        my $qual2 = $line[COUNT] > 1 ? getModalQual($line[QUAL2]) : $line[QUAL2];
        print "\@$name:2\n$line[SEQ2]\n+\n$qual2\n";  
    }
}
sub getModalQual {
    my %counts;    
    map { $counts{$_}++ } split("", $_[0]);
    my @quals = sort { $counts{$b} <=> $counts{$a} } keys %counts;
    $quals[0] x length($_[0]);
}
