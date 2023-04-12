use strict;
use warnings;

# extract read sequences for
#   a subset of single-alignment molecules for adapter model training
#   all SV molecules

# constants
use constant {
    QNAME => 0, 
    NODE1 => 1,
    NODE2 => 2,
    EDGE_TYPE => 3,
    _MAPQ => 4,
    SV_SIZE => 5,
    INSERT_SIZE => 6,
    XSTART => 7,
    XEND => 8,
    N_STRANDS => 9,
    #-------------
    ALIGNMENT => "A",
    JUNCTION  => "J"
};

# filter for qNames to retrieve
my ($nAlignments, $maxAlignments, $minSize, $nLines, $prevQName, $prevSize, %qNames) = (0, 10000, 500, 0);
while(my $line = <STDIN>){
    my @line = split("\t", $line);
    if($prevQName and $prevQName ne $line[QNAME]){
        if($nLines == 1){
            if($nAlignments <= $maxAlignments and $prevSize >= $minSize){
                $qNames{$prevQName} = ALIGNMENT;
                $nAlignments++;
            }
        } else {
            $qNames{$prevQName} = JUNCTION;
        }
        $nLines = 0;
    }
    $nLines++;
    $prevQName = $line[QNAME];
    $prevSize = $line[SV_SIZE];
}
$nLines > 1 and $qNames{$prevQName} = JUNCTION;

# run the fastq files once to retrieve the read sequences
my $fastqFiles = "$ENV{INPUT_DIR}/*.fastq.gz";
open my $inH,  "-|", "zcat $fastqFiles" or die "could not open: $fastqFiles\n";
while(my $qName = <$inH>){
    my $seq = <$inH>;
    my $discard = <$inH>;
    my $qual = <$inH>;
    $qName =~ m/^\@(\S+)/;
    $qNames{$1} and print join("\t", $qNames{$1}, $1, $seq); # $seq still has newline
}
close $inH;
