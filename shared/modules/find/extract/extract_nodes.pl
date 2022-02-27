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
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric sequence genome exclude);
resetCountFile();

# environment variables
fillEnvVar(\my $EXTRACT_PREFIX,  'EXTRACT_PREFIX');
fillEnvVar(\my $ACTION_DIR,      'ACTION_DIR');
fillEnvVar(\my $N_CPU,           'N_CPU');
fillEnvVar(\my $MIN_MAPQ,        'MIN_MAPQ');
fillEnvVar(\our $LIBRARY_TYPE,   'LIBRARY_TYPE');
fillEnvVar(\our $MIN_CLIP,       'MIN_CLIP');
fillEnvVar(\our $MAX_TLEN,       'MAX_TLEN');
fillEnvVar(\our $READ_LEN,       'READ_LEN');
use vars qw(%chromIndex);

# load additional dependencies
require "$ACTION_DIR/extract/parse_nodes.pl";

# initialize the genome
setCanonicalChroms();
initializeExclude($ENV{BAD_REGIONS_FILE}, $READ_LEN);

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
    _IS_PAIRED => 1,
    _PROPER => 2, # SAM FLAG bits
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
    #-------------
    MOL_ID => 0, # integer read-pair index
    MOL_CLASS => 1,
    MOL_STRAND => 2,
    IS_OUTER_CLIP1 => 3,
    IS_OUTER_CLIP2 => 4,
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    #-------------
    LEFT  => 0, # clip recovery direction for code readability
    RIGHT => 1,   
    #-------------
    IS_PROPER => 'P', # proper and anomalous molecule codes
    IS_SV => 'V',
    #-------------
    GAP        => 0, # SV evidence type codes, i.e. node classes
    SPLIT      => 1,
    OUTER_CLIP => 2,
};

# process data by molecule over multiple parallel threads
launchChildThreads(\&parseReadPairs);
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

# child process to parse BWA remap output
sub parseReadPairs {
    my ($childN) = @_;
    
    # working variables
    our ($excludeMol, $minMapQ, $alnN, $jxnN, @alns, @mol) = (0, 1e9, 1, 0);

    # output file handles
    my $spansFile = "$EXTRACT_PREFIX.spans.$childN.gz";
    open our $spansH, "|-", "gzip -c > $spansFile" or die "$error: could not open $spansFile: $!\n";
    my $nodesFile = "$EXTRACT_PREFIX.nodes.$childN.gz";
    open our $nodesH, "|-", "gzip -c > $nodesFile" or die "$error: could not open $nodesFile: $!\n";
  
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
                $mol[MOL_ID] = $aln[1];

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
