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
fillEnvVar(\our $EXTRACT_PREFIX,   'EXTRACT_PREFIX');
fillEnvVar(\our $ACTION_DIR,       'ACTION_DIR');
fillEnvVar(\our $READ_LEN,         'READ_LEN');
fillEnvVar(\our $MAX_INSERT_SIZE,  'MAX_INSERT_SIZE');
fillEnvVar(\our $MIN_SV_SIZE,      'MIN_SV_SIZE');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
map { require "$ACTION_DIR/extract/$_.pl" } qw(parse_nodes);
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svAmplicon";
map { require "$perlUtilDir/$_.pl" } qw(amplicons);
$perlUtilDir = "$ENV{MODULES_DIR}/parse_nodes";
map { require "$perlUtilDir/$_.pl" } qw(parse_nodes_support);
use vars qw(@amplicons);

# constants
use constant {
    END_READ_PAIR => '_ERP_',
    #-------------
    AMP_AMPLICON_ID => 0, # amplicon fields
    AMP_PROPER => 1,
    AMP_MOL_COUNT => 2,
    AMP_CHROM1 => 3,
    AMP_SIDE1 => 4,
    AMP_POS1 => 5,
    AMP_REF1 => 6,
    AMP_PRIMER1 => 7,
    AMP_CHROM2 => 8,
    AMP_SIDE2 => 9,
    AMP_POS2 => 10,
    AMP_REF2 => 11,
    AMP_PRIMER2 => 12,
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
    MERGE_FAILURE => "M",
    PROPER        => "P"
};

# working variables
our (@alns, @mol, $amplicon,
     @nodes,    @mapQs,    @cigars,    @alnQs,    @types,    @sizes,    @insSizes,    @outAlns,
     @alnNodes, @alnMapQs, @alnCigars, @alnAlnQs, @alnTypes, @alnSizes, @alnInsSizes, @alnAlns) = ();
my $anyUnmapped = 0;
our $APPEND_JXN_BASES = 1;

# process data by read pair
my ($prevPairName);
$| = 1;
while(my $line = <STDIN>){
    $nInputAlns++;
    chomp $line;
    my @aln = (split("\t", $line, 12))[QNAME..QUAL];     
    # name = molId:ampliconId:nOverlapBases:molCount:merged:readN 
    $aln[QNAME] =~ m/(.+):\d/; # strip the trailing readN to group read by readPair
    my $pairName = $1;
    if($prevPairName and $pairName ne $prevPairName){  
        parseReadPair();  
    }
    $aln[RNAME_INDEX] = $chromIndex{$aln[RNAME]} || 0; # non-canonical chroms are excluded, unmapped continue on (for now)   
    ($aln[FLAG] & _UNMAPPED) and $anyUnmapped = 1;
    @mol = split(":", $aln[QNAME]);    
    push @{$alns[$mol[READ_N] - 1]}, \@aln;
    $prevPairName = $pairName;
} 
parseReadPair(); 

# print summary information
printCount($nInputAlns, 'nInputAlns', 'input aligned segments over all read pairs');
printCount($nReadPairs, 'nReadPairs', 'input read pairs');

# parse bam read pairs
sub parseReadPair {
    $nReadPairs++;

    # reject molecules with even one unmapped read, can't be considered an amplicon match
    unless($anyUnmapped){
        $amplicon = $amplicons[$mol[AMPLICON_ID]];
        parseReadAlignments(READ1, @{$alns[READ1]});
        if($alns[READ2]){
            push @mapQs,    0;
            push @cigars,   "NA";
            push @alnQs,    0;
            push @types,    MERGE_FAILURE;
            push @sizes,    "NA";
            push @insSizes, "NA\tNA";
            push @outAlns,  [];
            parseReadAlignments(READ2, @{$alns[READ2]});
        }  
        fillJxnQs();
        printMolecule($mol[MOL_ID], @mol[AMPLICON_ID, MERGE_LEVEL, N_OVERLAP_BASES, IS_REFERENCE, MOL_COUNT]);
    }

    # reset for the next read-pair
    (@alns, @nodes, @mapQs, @cigars, @alnQs, @types, @sizes, @insSizes, @outAlns) = ();
    $anyUnmapped = 0;
}
