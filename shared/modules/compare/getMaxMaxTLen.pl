use strict;
use warnings;
use JSON;

# get a single, largest value for MAX_TLEN for making inter-sample SV comparisons

# initialize reporting
our $script = "maxMaxTLen";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
resetCountFile();

# environment variables
fillEnvVar(\my $FIND_GLOB_PREFIX, 'FIND_GLOB_PREFIX');
fillEnvVar(\my $TASK_DIR,         'TASK_DIR');
fillEnvVar(\my $GENOME,           'GENOME');

# parse sample-level information
my $maxMaxTLen = 0;
foreach my $sampleFile (glob("$FIND_GLOB_PREFIX.structural_variants.gz")){
    $sampleFile =~ m|$TASK_DIR/(.+)/|;
    my $sample = $1;   
    my $envFile = "$TASK_DIR/$sample/$sample.$GENOME.find.environment.json";
    my $env = decode_json( slurpFile($envFile) );
    my $maxTLen = $$env{MAX_TLEN}[0];
    printCount($maxTLen, "$sample-maxTLen", "$sample MAX_TLEN");
    $maxMaxTLen >= $maxTLen or $maxMaxTLen = $maxTLen;
}

# report the largest insert size in all libraries
print $maxMaxTLen;

1;

# # constants 
# use constant {
#     SV_ID => 1,
#     JUNCTION_NAME => 2,
#     MATCHING_NAMES => 3,
#     #-------------
#     TARGET_CLASS => 4,
#     JXN_TYPE => 5,
#     #-------------
#     CHROM_1 => 6,
#     SIDE_1 => 7,
#     POS_1 => 8,
#     CHROM_2 => 9,
#     SIDE_2 => 10,
#     POS_2 => 11,
#     #-------------
#     JXN_SEQ => 12,
#     MERGE_LEN => 13,
#     FAIDX_PADDING => 14,
#     GEN_REF_1 => 15,
#     GEN_REF_2 => 16,
#     #-------------
#     MICROHOM_LEN => 17,
#     MICROHOM_MATCH => 18,
#     JXN_BASES => 19,
#     SV_SIZE => 20 ,
#     #-------------
#     N_TOTAL => 21,
#     N_SPLITS => 22,
#     N_GAPS => 23,
#     N_OUTER_CLIPS => 24,
#     N_DUPLEX => 25,
#     N_DUPLEX_ALL => 26,
#     NET_STRAND_COUNT => 27,
#     NET_STRAND_COUNT_ALL => 28,
#     N_SHARED_PROPER => 29,
#     N_SHARED_PROPER_ALL => 30,
#     #------------- 
#     IS_MERGED => 31,
#     SEQ_LEN => 32,
#     UMI => 33,
#     MAPQ => 34,
#     #-------------
#     CHUNK_OFFSET => 35,
#     CHUNK_SIZE => 36,
#     #-------------
#     SAMPLE => 37,
#     # #-------------
#     # SV_ID => 0,
#     # MOL_IDS => 1, # comma-delimited
#     # #-------------
#     # TARGET_CLASS => 2, # TT tt etc.
#     # JXN_TYPE => 3,  # DIPTL
#     # FLAG1 => 4,
#     # RNAME1 => 5,
#     # INN_SIDE1 => 6,
#     # FLAG2 => 7,
#     # RNAME2 => 8,
#     # INN_SIDE2 => 9,
#     # #-------------
#     # PROX_OUT_POS1 => 10, # genome data 1
#     # PROX_JXN_POS1 => 11,
#     # PROX_JXN_POS_PLUS1 => 12,
#     # IS_CLIPPED1 => 13,
#     # JXN_JXN_POS1 => 14,
#     # JXN_MAPQ1 => 15,
#     # #-------------
#     # PROX_OUT_POS2 => 16, # genome data 2
#     # PROX_JXN_POS2 => 17,
#     # PROX_JXN_POS_PLUS2 => 18,
#     # IS_CLIPPED2 => 19,
#     # JXN_JXN_POS2 => 20,
#     # JXN_MAPQ2 => 21,
#     # #-------------
#     # CALL_LEN => 22, # junction call data
#     # JXN_SEQ => 23,
#     # JXN_DEPTH => 24,
#     # MICROHOM_LEN => 25,
#     # MICROHOM_MATCH => 26,
#     # JXN_BASES => 27,
#     # SV_SIZE => 28,
#     # #-------------
#     # N_TOTAL => 29, # junction evidence tallies
#     # N_SPLITS => 30,
#     # N_GAPS => 31,
#     # N_OUTER_CLIPS => 32,
#     # N_DUPLEX => 33,
#     # NET_STRAND_COUNT => 34,
#     # N_SHARED_PROPER => 35,
#     # N_UMI_PURGED => 36,
#     #-------------
#     SV_FILE => 37, # appended by compare.sh to identify the sample source
#     SAMPLE => 38, # added by this script based on SV_FILE
#     WAS_SEQUENCED => 39,
#     MAX_DIST => 40,
#     #=====================
#     LEFTWARD  => "L", # orientation of aligned source molecule relative to a mapped endpoint
#     RIGHTWARD => "R",
#     #---------------
#     INITIALIZED => 0,
#     PENDING => 1,
#     COMMITTED => 2
# };
