use strict;
use warnings;

# extract SV nodes from name-sorted read-pairs
# edges between nodes assembled later, in compile step

# initialize reporting
our $script = "extract_nodes";
our $error  = "$script error";
my ($nInputAlns, $nReadPairs, 
    $nThreadPairs, $nExcluded, $nFailedQuality) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms exclude);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\my $IS_SV_CAPTURE,     'IS_SV_CAPTURE', 1, 0);
fillEnvVar(\my $IS_USER_BAM,       'IS_USER_BAM');
fillEnvVar(\my $EXTRACT_PREFIX,    'EXTRACT_PREFIX');
fillEnvVar(\my $ACTION_DIR,        'ACTION_DIR');
fillEnvVar(\my $N_CPU,             'N_CPU');
fillEnvVar(\my $MIN_MAPQ,          'MIN_MAPQ');
fillEnvVar(\our $LIBRARY_TYPE,     'LIBRARY_TYPE');
fillEnvVar(\our $LIBRARY_ORIENTATION, 'LIBRARY_ORIENTATION', 1, "FR");
fillEnvVar(\our $MIN_CLIP,         'MIN_CLIP');
fillEnvVar(\our $MAX_TLEN,         'MAX_TLEN');
fillEnvVar(\our $READ_LEN,         'READ_LEN');
fillEnvVar(\our $BAD_REGIONS_FILE, 'BAD_REGIONS_FILE');

# load additional dependencies
require "$ACTION_DIR/extract/parse_nodes.pl";

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();
initializeExclude($BAD_REGIONS_FILE, $READ_LEN);

# constants
use constant {
    END_READ_PAIR => '_ERP_',
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
    QUAL => 10,
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
    IS_ORPHAN => 11,   # is this needed???
    INSERT_SIZE => 12, # is this needed???
    TARGET_CLASS => 13, # values added (or initialized) for svCapture
    SHARED_PROPER => 14, 
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    #-------------
    LEFT  => 0, # clip recovery direction for code readability
    RIGHT => 1,   
};

# working variables
my $isFF = ($LIBRARY_ORIENTATION eq "FF"); # else FR paired-end read orientation

# process data by molecule over multiple parallel threads
launchChildThreads(\&parseReadPair);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
while(my $line = <STDIN>){
    $nInputAlns++;
    my ($qName) = split("\t", $line, 2);    
    if($threadName and $qName ne $threadName){
        $nReadPairs++;        
        print $writeH END_READ_PAIR, "\t$nReadPairs\n";        
        $writeH = $writeH[$nReadPairs % $N_CPU + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $qName;
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
    
    # auto-flush output to prevent buffering and ensure proper feed to sort
    $| = 1;

    # working variables
    our ($excludeMol, $minMapQ, $alnN, $jxnN, @alns, @mol) = (0, 1e9, 1, 0);

    # output file handles
    my $spansFile = "$EXTRACT_PREFIX.spans.$childN.gz";
    open our $spansH, "|-", "gzip -c > $spansFile" or die "$error: could not open $spansFile: $!\n";
    my $nodesFile = "$EXTRACT_PREFIX.nodes.$childN.gz";
    open our $nodesH, "|-", "gzip -c > $nodesFile" or die "$error: could not open $nodesFile: $!\n";

    # # output file handles
    # my $endpointFile = "$ENV{EXTRACT_PREFIX}.endpoints.$childN.gz";
    # open our $endpointH, "|-", "gzip -c > $endpointFile" or die "$error: could not open $endpointFile: $!\n";
    # my $tLenFile = "$ENV{EXTRACT_PREFIX}.insertSizes.$childN.gz";
    # open our $tLenH, "|-", "gzip -c > $tLenFile" or die "$error: could not open $tLenFile: $!\n";
    # my $nodeFile = "$ENV{EXTRACT_PREFIX}.nodes.$childN.gz";
    # open our $nodeH, "|-", "gzip -c > $nodeFile" or die "$error: could not open $nodeFile: $!\n";
  
    # run aligner output one alignment at a time
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        my @aln = ((split("\t", $line, 11))[QNAME..SEQ], $alnN);     

        # parse output one source molecule at a time
        if($aln[0] eq END_READ_PAIR){
            $nThreadPairs++;
            
            # proceed if ALL alignment segments were of sufficient MAPQ
            # skip the entire molecule if ANY alignment segments were in an excluded region
            if(!$excludeMol and $minMapQ >= $MIN_MAPQ){

                # initialize molecule-level data based on bam source type
                if($IS_SV_CAPTURE){ # takes precedence over USER_BAM since realignment bam file forced by pipeline
                    @mol = (split(":", $alns[0][QNAME]), # MOL_ID to STRAND_COUNT2
                            (0) x 8);                    # MOL_CLASS to SHARED_PROPER
                    # $mol[MOL_STRAND] = $isTruSeq ? 0 : ($mol[STRAND_COUNT1] ? 0 : 1);
                } elsif($IS_USER_BAM){
                    @mol = ($aln[1],   # MOL_ID
                            (1) x 2,   # UMI1, UMI2 dummy values
                            ($alns[0][FLAG] & _IS_PAIRED) ^ _IS_PAIRED, # IS_MERGED as (0,1)
                            (0) x 11); # IS_DUPLEX to SHARED_PROPER
                } else { # genomex-mdi-tools align
                    @mol = (split(":", $alns[0][QNAME]), # MOL_ID to IS_MERGED
                            (0) x 11);                   # IS_DUPLEX to SHARED_PROPER
                }

                # identify pairs as proper or SV-containing and act accordingly

                # if($mol[IS_MERGED]){ # merged reads, expect just one alignment; any supplemental = SV
                #     if(@alns == 1 and !($alns[0][FLAG] & _UNMAPPED)){ # unmapped v. rare since remapping; dicard them
                #         parseMergedProper();
                #     } elsif(@alns > 1){
                #         parseMergedSplit();
                #     }                
                # } else { # unmerged reads, expect 2 alignments flagged as proper
                #     if(@alns == 2 and ($alns[0][FLAG] & _PROPER)){
                #         parseUnmergedProper();                    
                #     } elsif(@alns == 2 and
                #             !($alns[0][FLAG] & _UNMAPPED) and !($alns[1][FLAG] & _UNMAPPED)) {
                #         if($alns[0][FLAG] & _SUPPLEMENTARY or $alns[1][FLAG] & _SUPPLEMENTARY){
                #             $mol[IS_ORPHAN] = 1;
                #             parseMergedSplit();
                #         } else {
                #             parseUnmergedHiddenJunction();
                #         }
                #     } elsif(@alns >= 2) {
                #         parseUnmergedSplit();   
                #     } elsif(@alns == 1){ # one read of a pair failed consensus; treat remaining alignment as a single merged read
                #         $mol[IS_ORPHAN] = 1;
                #         parseMergedProper();
                #     }
                # }   

                # identify pairs as proper or SV-containing and act accordingly
                if($alns[0][FLAG] & _IS_PAIRED){ # unmerged reads, expect 2 alignments flagged as proper
                    if(@alns == 2 and ($alns[0][FLAG] & _PROPER)){ 
                        my ($read1, $read2) = ($alns[READ1][FLAG] & _FIRST_IN_PAIR) ? (READ1, READ2) : (READ2, READ1);
                        commitProperMolecule($alns[$read1], $alns[$read2], LEFT,  RIGHT); 
                    } elsif(@alns == 2 and !($alns[0][FLAG] & _UNMAPPED) and !($alns[1][FLAG] & _UNMAPPED)) {
                        parseUnmergedHiddenJunction();
                    } else {
                        parseUnmergedSplit();   
                    }
                } else { # merged or orphaned reads, expect just one alignment; any supplemental = SV
                    if(@alns == 1){ # equivalent to a _PROPER flag for merged/singleton reads
                        $alns[MERGED_READ][FLAG] & _UNMAPPED or  # unmapped singleton reads are discarded
                            commitProperMolecule($alns[MERGED_READ], $alns[MERGED_READ], RIGHT, LEFT);
                    } else {
                        parseMergedSplit();
                    }                  
                } 
            } elsif($excludeMol) {
                $nExcluded++;
            } else {
                $nFailedQuality++;
            }
    
            # prepare for next read-pair
            $excludeMol = 0;
            $alnN = 1;
            $jxnN = 0;
            @alns = ();
            $minMapQ = 1e9;

        } elsif(!$excludeMol){ # add new alignment to growing source molecule
            $aln[RNAME] = $chromIndex{$aln[RNAME]} || 0; # non-canonical chroms are excluded, unmapped continue on (for now)   
            $excludeMol = isExcludedPosition(@aln[RNAME, POS]) and next;
            push @alns, \@aln;
            $minMapQ <= $aln[MAPQ] or $minMapQ = $aln[MAPQ];
            $alnN++; # unique identifier to establish alignment edges between nodes   
        }
    }

    # close outputs
    close $spansH;
    close $nodesH;

    # report a thread summary to log
    printCount($nThreadPairs,   "nThreadPairs-$childN",   "total read pairs processed in thread $childN");
    printCount($nExcluded,      "nExcluded-$childN",      "rejected read pairs in excluded regions in thread $childN");
    printCount($nFailedQuality, "nFailedQuality-$childN", "rejected read pairs failed MIN_MAPQ $MIN_MAPQ in thread $childN");
}

# UNMERGED READS
# 1111111^1111111>-------^----------------^---------
# -------^---------------^---------<222222^222222222

# MERGED READS
# -------^---------------^----------------^---------
# -------^---------------^----------------^---------
