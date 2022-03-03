use strict;
use warnings;

# extract count/distribution handling
# creates the files that described the source molecules, per thread

# constants
use constant {
    IS_PROPER => 'P', # proper and anomalous molecule classes
    IS_SV     => 'V',
    #-------------
    MAX_FAMILY_COUNT => 100
};

# operating parameters
use vars qw($error $IS_COLLATED $isTargeted $EXTRACT_PREFIX $MAX_TLEN);
our $maxInsertSize = int($MAX_TLEN * 1.25);  
our $isCountStrands = ($IS_COLLATED and $isTargeted);   
our (%insertSizes, %strandCounts);

# print proper source molecule insert sizes by targetClass
# used in all extractions runs
sub printInsertSizeFile {
    my ($childN) = @_;
    my $file = getFileStream('insert_sizes', $childN);
    open my $crosstabH, "|-", $$file{stream} or die $$file{error};
    my @insertSizes = 0..$maxInsertSize;
    my @targetClasses = ('X', 'TT', 'TA', 'AA', '--'); # the target classes that could support proper molecules
    print $crosstabH join("\t", "insertSize", 'untargeted', 'TT', 'TA', 'AA', 'xx'), "\n"; # R-compatible column names
    foreach my $insertSize(@insertSizes){
        print $crosstabH "$insertSize";
        foreach my $targetClass(@targetClasses){
            my $count = $insertSizes{$targetClass}[$insertSize] || 0;
            print $crosstabH "\t$count";
        }
        print $crosstabH "\n";
    }
    close $crosstabH;     
}

# print strand family sizes by targetClass and molClass
# only applicable if collated and targeted
sub printCrosstabFile {
    my ($childN) = @_;
    $isCountStrands or return;
    my $file = getFileStream('strand_counts', $childN);
    open my $crosstabH, "|-", $$file{stream} or die $$file{error};
    my @tableCounts = 0..MAX_FAMILY_COUNT;
    my @types;
    foreach my $targetClass ('TT', 'TA', 'AA', 'tt', 'ta', 'aa', 't-', 'a-', '--'){
        foreach my $molClass (IS_PROPER, IS_SV){ # not all combinations are possible, but that's fine
            push @types, "$targetClass\t$molClass";
        }
    }
    print $crosstabH join("\t", qw(targetClass molClass count1), @tableCounts), "\n";
    foreach my $type(@types){
        foreach my $count1(@tableCounts){
            print $crosstabH "$type\t$count1";
            foreach my $count2(@tableCounts){
                print $crosstabH "\t", ($strandCounts{$type} and $strandCounts{$type}[$count1]) ?
                                       ($strandCounts{$type}[$count1][$count2] || 0) : 0;
            }
            print $crosstabH "\n";
        }
    }
    close $crosstabH;     
}

1;
