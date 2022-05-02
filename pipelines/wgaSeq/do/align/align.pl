use strict;
use warnings;

#----------------------------------------------------------------
# align a series of individual cells to genome and write output to a single stream
#----------------------------------------------------------------

# initialize reporting
our $action = 'align';
my ($nCells) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
require "$perlUtilDir/workflow.pl";
resetCountFile();

# environment variables
fillEnvVar(\my $nCpu,       'N_CPU');
fillEnvVar(\my $totalRamInt,'TOTAL_RAM_INT');
fillEnvVar(\my $tmpDirWrk,  'TMP_DIR_WRK');
fillEnvVar(\my $logPrefix,  'LOG_FILE_PREFIX');
fillEnvVar(\my $actionDir,  'ACTION_DIR');
fillEnvVar(\my $genomeFasta,'BWA_GENOME_FASTA');
fillEnvVar(\my $cellIds,    'CELL_IDS');
fillEnvVar(\my $inputMode,  'INPUT_MODE');

# set memory usage
my $alignRam = 5e9; # RAM set aside for BWA and all other non-sorting programs in stream
my $sortRam = int(($totalRamInt - $alignRam) / 2); # distribute remaining RAM over two sort processes
my $sortRamPerCpu = int($sortRam / $nCpu);

# constants
use constant {
    QNAME => 0, # sam fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    QUAL => 10,
    #-------------
    _RNAME => 0, # grouping stream columns
    _GROUP_POS1 => 1,
    _GROUP_POS2 => 2,
    _POS => 3,
    _MAPQ => 4,
    _TLEN => 5,
    _DUP_RATE => 6,
    #-------------
    _PAIRED => 1,
    _PROPER_PAIR => 2,
    _PROPER_PAIR_BITS => 1 + 2,
};

# get inputs
my @cellIds = split(/\s+/, $cellIds);

# initialize outputs
my $bamListFile  = "$logPrefix.$action.bam.list";
my $cellListFile = "$logPrefix.$action.cell.list";
open my $bamListH, ">", $bamListFile  or die "could not open $bamListFile for writing\n";
open my $cellH,    ">", $cellListFile or die "could not open $cellListFile for writing\n";
print $cellH join(",", qw(barcode cell_id)), "\n";
my $bwaLogFile = "$logPrefix.$action.bwa.log";
unlink $bwaLogFile;

# collect a common bam header for all cell output files
sub getCellFastqs {
    my ($cellId) = @_; # working from within INPUT_DIR
    glob($inputMode eq "directory" ? "$cellId/*.fastq.gz" : $cellId."_*.fastq.gz");
}
my $fastQs = join(" ", getCellFastqs($cellIds[0]));
my $bwa = "bwa mem -Y -t $nCpu $genomeFasta $fastQs 2>/dev/null";
my $samtools = "samtools view -H";
my $bamHeader = qx/$bwa | $samtools/;

# align one cell at a time with parallelization
foreach my $cellId(@cellIds){
    my $time = localtime();
    
    # add cell to cell list in abbreviated 10x CellRanger format    
    $nCells++;
    my $cellN = $nCells - 1;
    print STDERR "  $cellId (cell #$cellN), started $time\n";
    print $cellH join(",", $cellId, $cellN), "\n";
    
    # set files
    my @fastQs = getCellFastqs($cellId);
    @fastQs == 2 or die "did not get expected two FASTQ files\n".join("\n", @fastQs)."\n";
    my $tmpBam = "$tmpDirWrk/cell_$cellN.sorted.bam";
    print $bamListH "$tmpBam\n";
    -e $tmpBam and next; # here to help with job recovery without remapping prior successes
    
    # align reads and remove duplicates
    # output does not retain full alignments, just the needed information about unique source molecules
    my $fastpLogPrefix = "$logPrefix.$cellId.fastp";
    my $fastp = "fastp --thread $nCpu --in1 $fastQs[0] --in2 $fastQs[1] --stdout ".
                "--merge --include_unmerged --disable_quality_filtering ".
                "--html $fastpLogPrefix.html --json $fastpLogPrefix.json --report_title \"$cellId\" 2>/dev/null";
    my $bwa      = "bwa mem -p -Y -t $nCpu $genomeFasta - 2>>$bwaLogFile";    
    my $parseBWA = "perl $actionDir/align/parse_bwa.pl";
    my $sort     = "sort --parallel=$nCpu -T $tmpDirWrk -S $sortRam"."b --compress-program=pigz -k1,1 -k2,2n -k3,3n";
    my $group    = "bedtools groupby -g 1,2,3 -c 4,5,6,6 -o first,max,first,count";
    my $inStream = "$fastp | $bwa | $parseBWA | $sort | $group";
    open my $inH, "-|", $inStream or die "could not open bwa stream: $!\n";    

    # process to bam for subsequent merging
    my $toBam = "samtools view --threads $nCpu -b -";
    my $bamSort = "samtools sort --threads $nCpu -m $sortRamPerCpu ".
                  "-T $tmpDirWrk/align_sort - 2>/dev/null"; # yes, we DO need to resort by POS (not groupPos)           
    my $outStream = "$toBam | $bamSort > $tmpBam";    
    open my $outH, "|-", $outStream or die "could not open samtools stream: $!\n";

    # process to sam with 10x scCNV-compatible CB:Z: tag
    my $molId = 1;
    print $outH $bamHeader;
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        print $outH join("\t",
            $molId, 
            _PROPER_PAIR_BITS,
            @f[_RNAME, _POS, _MAPQ],
            '1M', # samtools demands a CIGAR string, doesn't matter what it is, we will never use it
            '=',
            1,
            $f[_TLEN],
            '*',
            '*',
            "XC:i:$f[_DUP_RATE]",
            "CB:Z:$cellId"
        ), "\n";        
        $molId++;
    }
    
    # clean up
    close $inH;
    close $outH;
}
    
# report counts
printCount($nCells, 'nCells', 'number of cells aligned to genome');

# clean up
close $cellH;
close $bamListH;
