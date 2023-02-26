use strict;
use warnings;

# extract SV information from collated, name-sorted amplicon read-pairs

# initialize reporting
our $script = "extract_nodes";
our $error  = "$script error";
my ($nInputAlns, $nReadPairs) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $COLLATE_PREFIX,   'COLLATE_PREFIX');
fillEnvVar(\our $ACTION_DIR,       'ACTION_DIR');
fillEnvVar(\our $N_CPU,            'N_CPU'); # user options, or derived from them
fillEnvVar(\our $READ_LEN,         'READ_LEN');
fillEnvVar(\our $MAX_INSERT_SIZE,  'MAX_INSERT_SIZE');

# load additional dependencies
map { require "$ACTION_DIR/extract/$_.pl" } qw(parse_nodes);

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# constants
use constant {
    END_READ_PAIR => '_ERP_',
    #-------------
    AMP_CHROM1 => 0, # amplicon fields
    AMP_SIDE1 => 1,
    AMP_POS1 => 2,
    AMP_CHROM2 => 3, # amplicon fields
    AMP_SIDE2 => 4,
    AMP_POS2 => 5,
    AMP_MOL_COUNT => 6,
    AMP_AMPLICON_ID => 7,
    AMP_INDEX => 8,
    AMP_PROPER => 9,
    #-------------
    QNAME => 0, # SAM fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    ALN_N => 10,
    RNAME_INDEX => 11,
    #-------------
    _IS_PAIRED => 1, # SAM FLAG bits
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
    #-------------
    MOL_ID => 0,    # molecule-level data, carried here in QNAME
    AMPLICON_ID => 1, 
    N_OVERLAP_BASES => 2, 
    MOL_COUNT => 3, 
    IS_MERGED => 4, 
    READ_N => 5, 
    MOL_CLASS => 5, # values added by extract_nodes regardless of pipeline or bam source
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0
};

# load the pre-determined amplicons
our @amplicons;
open my $inH, "<", "$COLLATE_PREFIX.amplicons.txt" or die "missing file: $COLLATE_PREFIX.amplicons.txt\n";
while(my $line = <$inH>){
    chomp $line;
    my @line = split("\t", $line);
    $amplicons[$line[AMP_AMPLICON_ID]] = \@line;
}
close $inH;

# process data by molecule over multiple parallel threads
launchChildThreads(\&parseReadPair);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
while(my $line = <STDIN>){
    $nInputAlns++;
    my ($qName) = split("\t", $line, 2);  
    # name = molId:ampliconId:nOverlapBases:molCount:merged:readN 
    $qName =~ m/(.+):\d/; # strip the trailing readN to group read by readPair
    my $pairName = $1;
    if($threadName and $pairName ne $threadName){
        $nReadPairs++;        
        print $writeH END_READ_PAIR, "\t$nReadPairs\n";        
        $writeH = $writeH[$nReadPairs % $N_CPU + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $pairName;
}
$nReadPairs++;
print $writeH END_READ_PAIR, "\t$nReadPairs\n";      
finishChildThreads();

# print summary information
printCount($nInputAlns, 'nInputAlns', 'input aligned segments over all read pairs');
printCount($nReadPairs, 'nReadPairs', 'input read pairs');

# child process to parse bam read pairs
sub parseReadPair {
    my ($childN) = @_;

    # working variables
    our ($alnN, $jxnN, @alns, @mol, $amplicon) = (1, 0);
  
    # run aligner output one alignment at a time
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        my @aln = ((split("\t", $line, 11))[QNAME..SEQ], $alnN);     

        # parse output one source molecule at a time
        if($aln[0] eq END_READ_PAIR){

            # reject molecules with even one unmapped read
            my $anyUnmapped = 0;
            foreach my $aln (@alns){ ($$aln[FLAG] & _UNMAPPED) and $anyUnmapped = 1 } 
            if(!$anyUnmapped){

                # initialize molecule-level data 
                @mol = split(":", $alns[0][QNAME]); # has READ_N in eventual MOL_CLASS slot
                $amplicon = $amplicons[$mol[AMPLICON_ID]];

                # identify pairs as proper or SV-containing and act accordingly
                if($mol[IS_MERGED]){ # merged reads, expect just one alignment; any supplemental = SV
                    if(@alns == 1){
                        $alns[MERGED_READ][FLAG] & _UNMAPPED or # unmapped singleton reads are discarded
                            commitContiguousAlignment($alns[MERGED_READ]);
                    } else {
                        parseMergedSplit();
                    } 
                } else { # unmerged read pairs, expect 2 alignments
                    foreach my $aln (@alns){ # set the SAM FLAG pairing bits for unmerged reads; PROPER flag not meaningful and not used
                        my @read = split(":", $$aln[QNAME]);
                        $$aln[FLAG] |= _IS_PAIRED + ($read[READ_N] == 1 ? _FIRST_IN_PAIR : _SECOND_IN_PAIR);
                    }
                    if(@alns == 2) {
                        parseUnmergedHiddenJunction();
                    } else {
                        parseUnmergedSplit();   
                    }  
                }                 
            }

            # prepare for next read-pair
            $alnN = 1;
            $jxnN = 0;
            @alns = ();

        } else{ # add new alignment to growing source molecule
            $aln[RNAME_INDEX] = $chromIndex{$aln[RNAME]} || 0; # non-canonical chroms are excluded, unmapped continue on (for now)   
            push @alns, \@aln;
            $alnN++; # unique identifier to establish alignment edges between nodes   
        }
    }
}
