use strict;
use warnings;

# extract comprehensive CNV information from name-sorted bam stream
# for a given, single cell

# initialize reporting
our $script = "extract";
our $error  = "$script error";
my ($nReadPairs, $nDiscarded) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $ALIGNMENT_DIR,       'ALIGNMENT_DIR');
fillEnvVar(\our $ALIGNMENT_FILE_TYPE, 'ALIGNMENT_FILE_TYPE');
fillEnvVar(\our $DATA_NAME,     'DATA_NAME');
fillEnvVar(\our $GENOME,        'GENOME');
fillEnvVar(\our $ACTION_DIR,    'ACTION_DIR');
fillEnvVar(\our $N_CPU,         'N_CPU');
fillEnvVar(\our $MIN_MAPQ,      'MIN_MAPQ');
fillEnvVar(\our $MODULES_DIR,   'MODULES_DIR');

# parse target cell
my ($alignmentFile) = @ARGV;
my $cellId = $alignmentFile;
$cellId =~ s/^$ALIGNMENT_DIR\///;
$cellId =~ s/^$DATA_NAME\.//;
$cellId =~ s/^$GENOME\.//;
$cellId =~ s/\.$ALIGNMENT_FILE_TYPE$//;
$cellId =~ s/\.name$//;
our $READ_LEN = qx(samtools view -f 1 $alignmentFile | head -n 1 | awk '{print length(\$10)}');
chomp $READ_LEN;
our $MAX_TLEN = 2000;
# our $MAX_TLEN = qx(samtools view $alignmentFile | head -n 100000 | perl $MODULES_DIR/library/get_max_TLEN.pl);
# chomp $MAX_TLEN;

# load additional cell-dependent dependencies
map { require "$ACTION_DIR/$_.pl" } qw(parse_nodes);

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# constants
use constant {
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNAME_INDEX => 6,
    #-------------
    _IS_PAIRED => 1, # SAM FLAG bits
    _PROPER => 2,
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
    #-------------
    MOL_ID => 0,    # molecule-level data, carried here in QNAME
    UMI1 => 1,      # appended to QNAME by genomex-mdi-tools align / svCapture make_consensus.pl
    UMI2 => 2,      #     will be added by extract_nodes if user-bam was prepared elsewhere
    IS_MERGED => 3,
    IS_DUPLEX => 4, # appended to QNAME by svCapture make_consensus.pl (but not genomex-mdi-tools align)
    STRAND_COUNT1 => 5,
    STRAND_COUNT2 => 6,
    MOL_CLASS => 7, # values added by extract_nodes regardless of pipeline or bam source
    MOL_STRAND => 8,
    IS_OUTER_CLIP1 => 9,
    IS_OUTER_CLIP2 => 10,
    TARGET_CLASS => 11, # values added (or initialized) for svCapture only
    SHARED_PROPER => 12, 
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    #-------------
    LEFT  => 0, # clip recovery direction for code readability
    RIGHT => 1,   
    #-------------
    IS_PROPER => 'P', # proper and anomalous molecule classes
    IS_SV     => 'V',
    #-------------
    BIN_SIZE => 20000 # fixed to match 10x scCNV
};

# parse aligned reads one source molecule at a time
our (@alns, @jxns, $molKey);
my ($minMapQ, $sumChromIndex, $prevQName) = (1e9, 0);
open my $inH, "-|", "slurp -s 250M $alignmentFile | samtools view -" or 
    die "could not open: $alignmentFile\n";
while(my $line = <$inH>){
    chomp $line;
    my @aln = (split("\t", $line, 7))[QNAME..CIGAR]; 
    if($prevQName and $prevQName ne $aln[QNAME]){
        parseReadPair();
        (@alns, @jxns, $molKey) = ();        
        ($minMapQ, $sumChromIndex) = (1e9, 0);
    } 
    $aln[RNAME_INDEX] = $chromIndex{$aln[RNAME]} || 0; # non-canonical chroms are excluded, unmapped continue on (for now) 
    $sumChromIndex += $aln[RNAME_INDEX];
    $minMapQ <= $aln[MAPQ] or $minMapQ = $aln[MAPQ];
    push @alns, \@aln;
    $prevQName = $aln[QNAME];
}
parseReadPair();
close $inH;

# print summary information
my $keptPairs = commify($nReadPairs - $nDiscarded);
$nReadPairs = commify($nReadPairs);
print STDERR "$keptPairs of $nReadPairs input read pairs kept for cell $cellId\n";

# child process to parse bam read pairs
sub getCoverageBin {
    my ($chromI, $pos) = @_;
    join(":", $chromI, int(($pos - 1) / BIN_SIZE) + 1);
}
sub parseReadPair {
    $nReadPairs++;

    # proceed if ALL alignment segments were of sufficient MAPQ
    # and at least one alignment was on a canonical chromosome
    if($minMapQ >= $MIN_MAPQ and $sumChromIndex > 0){

        # flip the called strand of all alignments arising from the 2nd read of an FR pair (even orphans)
        # from this point forward, all libraries are handled as FF libraries
        foreach my $aln(@alns){
            if($$aln[FLAG] & _SECOND_IN_PAIR){
                $$aln[FLAG] ^= _REVERSE;
            }
        }

        # identify pairs as proper or SV-containing and act accordingly
        if($alns[0][FLAG] & _IS_PAIRED){ # unmerged read pairs, expect 2 alignments flagged as proper
            if(@alns == 2 and ($alns[0][FLAG] & _PROPER)){ 
                my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
                commitProperMolecule($alns[$read1], $alns[$read2]); 
            } elsif(@alns == 2 and !($alns[0][FLAG] & _UNMAPPED) and !($alns[1][FLAG] & _UNMAPPED)) {
                parseUnmergedHiddenJunction();
            } else {
                parseUnmergedSplit();   
            }
        } else { # merged or orphaned reads, expect just one alignment; any supplemental = SV
            if(@alns == 1){ # equivalent to a _PROPER flag for merged/singleton reads
                $alns[MERGED_READ][FLAG] & _UNMAPPED or # unmapped singleton reads are discarded
                    commitProperMolecule($alns[MERGED_READ], $alns[MERGED_READ]);
            } else {
                parseMergedSplit();
            }                  
        } 

        # commit bins and junction edges in a single row with a molecule identifier based on outer endpoints
        if(@alns){
            my %bins;
            foreach my $aln(@alns){
                $bins{getCoverageBin($$aln[RNAME_INDEX], $$aln[POS])} = 1;
                $bins{getCoverageBin($$aln[RNAME_INDEX], getEnd($$aln[POS], $$aln[CIGAR]))} = 1;
            }
            print join("\t", $molKey, scalar(@alns), join("::", keys %bins), join(":::", @jxns) || "X"), "\n";           
        }

        #########################
        # $nReadPairs >= 10000 and die "----------------\n"; 

    } else {
        $nDiscarded++;
    }
}
