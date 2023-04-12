use strict;
use warnings;

# extract junction information from name-sorted amplicon long reads

# initialize reporting
our $script = "extract_nodes";
our $error  = "$script error";
my ($nInputAlns, $nMolecules) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $EXTRACT_PREFIX,   'EXTRACT_PREFIX');
fillEnvVar(\our $ACTION_DIR,       'ACTION_DIR');
our $N_CPU = 1; # at present, large write blocks create problems with auto-flush to output stream
our $MIN_SV_SIZE = 1e6; # thereby disallowing all junctions within individual CIGAR strings

# load additional dependencies
map { require "$ACTION_DIR/extract/$_.pl" } qw(parse_nodes);
$perlUtilDir = "$ENV{MODULES_DIR}/parse_nodes";
map { require "$perlUtilDir/$_.pl" } qw(parse_nodes_support);

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# constants
use constant {
    END_MOLECULE => '_ENDMOL_',
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
    MERGE_LEVEL => 2,
    N_OVERLAP_BASES => 3, 
    IS_REFERENCE => 4,
    MOL_COUNT => 5, 
    READ_N => 6, 
    # MOL_CLASS => 6, # values added by extract_nodes regardless of pipeline or bam source
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    #-------------
    LEFT  => 0, # clip recovery direction for code readability
    RIGHT => 1,  
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    PROPER        => "P"
};

# process data by molecule over multiple parallel threads
launchChildThreads(\&parseReadPair);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
while(my $line = <STDIN>){
    $nInputAlns++;
    my ($qName) = split("\t", $line, 2);  
    if($threadName and $qName ne $threadName){
        $nMolecules++;        
        print $writeH END_MOLECULE, "\t$nMolecules\n";        
        $writeH = $writeH[$nMolecules % $N_CPU + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $qName;
}
$nMolecules++;
print $writeH END_MOLECULE, "\t$nMolecules\n";      
finishChildThreads();

# print summary information
printCount($nInputAlns, 'nInputAlns', 'input aligned segments over all molecules');
printCount($nMolecules, 'nMolecules', 'input molecules');

# child process to parse bam read pairs
sub parseReadPair {
    my ($childN) = @_;
    
    # auto-flush output to prevent buffering and ensure proper feed
    $| = 1;

    # working variables
    our (@alns, 
         @nodes,    @types,    @mapQs,    @sizes,    @insSizes,    @outAlns,
         @alnNodes, @alnTypes, @alnMapQs, @alnSizes, @alnInsSizes, @alnAlns) = ();
    my $anyUnmapped = 0;
  
    # run aligner output one alignment at a time
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        my @aln = (split("\t", $line, 12))[QNAME..QUAL];     

        # parse output one source molecule at a time
        if($aln[0] eq END_MOLECULE){

            # reject molecules with even one unmapped read, can't be considered and amplicon match
            unless($anyUnmapped){
                parseReadAlignments(@alns); 
                fillJxnMapQs();
                printMolecule($aln[QNAME]);
            }

            # prepare for next read-pair
            (@alns, @nodes, @types, @mapQs, @sizes, @insSizes, @outAlns) = ();
            $anyUnmapped = 0;

        } else{ # add new alignment to growing source molecule
            $aln[RNAME_INDEX] = $chromIndex{$aln[RNAME]} || 0; # non-canonical chroms are excluded, unmapped continue on (for now)   
            ($aln[FLAG] & _UNMAPPED) and $anyUnmapped = 1;
            push @alns, \@aln;
        }
    }
}
