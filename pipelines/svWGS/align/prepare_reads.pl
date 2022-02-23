use strict;
use warnings;

# action:
#     prepare interleaved FASTQ from different input types, including SRA files
# expects:
#     source $MODULES_DIR/scan/set_read_file_vars.sh (sets FASTQ_FILE1, FASTQ_FILE2, SRA_FILES)
#     input as either paired fastq.gz files or a set of .sra files

# initialize reporting
our $action  = "prepare_reads";
my $nReadPairs = 0;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
require "$perlUtilDir/workflow.pl";
resetCountFile();

# environment variables
fillEnvVar(\my $fastq1,  'FASTQ_FILE1');
fillEnvVar(\my $fastq2,  'FASTQ_FILE2');
fillEnvVar(\my $sraFiles,'SRA_FILES');

# constants
use constant {
    QNAME => 0,
    SEQ   => 1,
    QUAL  => 2
};

# set the file input handles
if($fastq1){
    open my $inH1, "-|", "zcat $fastq1" or throwError("could not open $fastq1: $!");
    open my $inH2, "-|", "zcat $fastq2" or throwError("could not open $fastq2: $!");
    runReadPairs($inH1, $inH2);
    close $inH1;
    close $inH2;    
} else {
    foreach my $sraFile(split(" ", $sraFiles)){
        open my $inH, "-|", "fastq-dump --stdout --split-files $sraFile" or throwError("could not open $sraFile: $!");
        runReadPairs($inH, $inH);
        close $inH;        
    }
}

# run the paired reads
# no need to multi-thread since the downstream alignment process is rate limiting
sub runReadPairs {
    my ($inH1, $inH2) = @_;
    while (my $read1 = getRead($inH1)){
        my $read2 = getRead($inH2);
        $nReadPairs++;
        print join("\n", $$read1[QNAME], $$read1[SEQ], '+', $$read1[QUAL]), "\n"; # print interleaved read pairs
        print join("\n", $$read1[QNAME], $$read2[SEQ], '+', $$read2[QUAL]), "\n";        
    }    
}

# parse a FASTQ set of 4 lines
sub getRead {
    my ($inH) = @_;

    # name line
    my $name = <$inH>;
    $name or return;
    chomp $name;
    my @name = split(/\s/, $name, 2); 
    
    # seq line  
    my $seq = <$inH>;
    chomp $seq;
    
    # discard line
    my $discard = <$inH>;

    # qual line
    my $qual = <$inH>;
    chomp $qual;
    
    return [$name[0], $seq, $qual];
}

# print summary information
printCount($nReadPairs,   'nReadPairs',   'total input read pairs');
