use strict;
use warnings;

# extract read sequences and base qualities for
#   a subset of single-alignment molecules for adapter model training
#   all SV molecules

# constants
use constant {
    QNAME => 0,
    # NODE1 => 1,
    # CIGAR1 => 2,
    # BLAST_IDENTITY1 => 3,
    # GAP_COMPRESSED_IDENTITY1 => 4,
    # NODE2 => 5,
    # CIGAR2 => 6,
    # BLAST_IDENTITY2 => 7,
    # GAP_COMPRESSED_IDENTITY2 => 8,
    # EDGE_TYPE => 9,
    # _MAPQ => 10,
    # SV_SIZE => 11,
    # INSERT_SIZE => 12,
    # XSTART => 13,
    # XEND => 14,
    # EDGE_CLASS => 15,
    # N_STRANDS => 16,
    #-------------
    ALIGNMENT => "A",
    JUNCTION  => "J"
};

# filter for qNames to retrieve
my ($nLines, $prevQName, %qNames) = (0);
while(my $line = <STDIN>){
    my @line = split("\t", $line, 2);
    if($prevQName and $prevQName ne $line[QNAME]){
        $qNames{$prevQName} = $nLines == 1 ? ALIGNMENT : JUNCTION;
        $nLines = 0;
    }
    $nLines++;
    $prevQName = $line[QNAME];
}
$qNames{$prevQName} = $nLines == 1 ? ALIGNMENT : JUNCTION;

# run the fastq files once to retrieve the required read sequences and quals
my $fastqFiles = "$ENV{INPUT_DIR}/*.fastq.gz";
open my $inH,  "-|", "zcat $fastqFiles" or die "could not open: $fastqFiles\n";
open my $idxH, ">", "$ENV{SEQUENCES_FILE}.index" or die "could not open file: $ENV{SEQUENCES_FILE}.index\n";
my $offset = 0;
while(my $qName = <$inH>){
    my $seq = <$inH>;
    my $discard = <$inH>;
    my $qual = <$inH>;    
    $qName =~ m/^\@(\S+)/;
    if($qNames{$1}){
        # chomp $seq;
        # my $line = join("\t", $qNames{$1}, $1, $seq, $qual);  # $qual still has newline

        my $line = join("\t", $qNames{$1}, $1, $seq);  # $seq still has newline

        print $line;
        my $nChar = length($line);   
        print $idxH join("\t", $1, $offset, $nChar), "\n";
        $offset += $nChar;        
    }
}
close $inH;
close $idxH;
