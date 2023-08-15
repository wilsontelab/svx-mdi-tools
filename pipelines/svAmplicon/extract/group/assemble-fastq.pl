use strict;
use warnings;

# write a FASTQ stream of kept and grouped molecule sequences suitable for re-alignment
# each molecule may have one merged or two unmerged reads
# name = molId:ampliconId:merged:nOverlapBases:isReference:molCount:readN

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(numeric);

use constant {
    AMPLICON=> 0, # the grouped sequences shared across one or more molecules
    SEQ1    => 1,
    SEQ2    => 2,
    #-------------
    QUAL1   => 3, # a comma-separated concatenation of all individual read qualities
    QUAL2   => 4,
    MERGED  => 5,
    OVERLAP => 6,
    MOL_ID  => 7, # one representative molecule of a given sequence
    IS_REFERENCE => 8,
    COUNT   => 9  # number of times this specific sequence was encountered, with all any N bases, sequencing errors, etc.
};

# parse the stream
while(my $line = <STDIN>){ # each line summarizes one unique encountered read sequence
    chomp $line;
    my @line = split("\t", $line);
    my $name = join(":", @line[MOL_ID, AMPLICON, MERGED, OVERLAP, IS_REFERENCE, COUNT]);
    my $qual1 = $line[COUNT] > 1 ? getMaxQual($line[QUAL1]) : $line[QUAL1];
    print "\@$name:1\n$line[SEQ1]\n+\n$qual1\n";
    if($line[SEQ2] ne "*"){
        my $qual2 = $line[COUNT] > 1 ? getMaxQual($line[QUAL2]) : $line[QUAL2];
        print "\@$name:2\n$line[SEQ2]\n+\n$qual2\n";  
    }
}
sub getMaxQual { # for multiply detected reads, return the maximum observed base quality at each position over all reads
    my ($readQuals) = @_;
    my @readQuals = split(",", $readQuals);
    my @maxBaseQuals;
    foreach my $baseN(1..length($readQuals[0])){ # perfectly fine to do alphanumeric sort on phred characters
        my @baseQuals = sort { $b cmp $a } map { substr($readQuals[$_], $baseN - 1, 1) } 0..$#readQuals;
        push @maxBaseQuals, $baseQuals[0];
    }
    join("", @maxBaseQuals);
    # my %counts;    # TODO: use better aggregate QUAL, e.g., via $val = ord($base) - 33 ??
    # map { $counts{$_}++ } split("", $_[0]);
    # my @quals = sort { $counts{$b} <=> $counts{$a} } keys %counts;
    # $quals[0] x length($_[0]);
}
