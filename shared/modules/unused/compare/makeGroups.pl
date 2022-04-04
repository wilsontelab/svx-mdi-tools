use strict;
use warnings;
use JSON;

# break SVs into smaller continuity groups prior to finding collision
# also add an SV key (sample + SV ID) for tracking comparisons

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);

# environment variables
fillEnvVar(\my $MAX_MAX_TLEN, 'MAX_MAX_TLEN');

# constants 
use constant {
    SV_ID => 0,
    JUNCTION_NAME => 1,
    MATCHING_NAMES => 2,
    #-------------
    TARGET_CLASS => 3,
    JXN_TYPE => 4,
    #-------------
    CHROM_1 => 5,
    SIDE_1 => 6,
    POS_1 => 7,
    CHROM_2 => 8,
    SIDE_2 => 9,
    POS_2 => 10,
    #-------------
    JXN_SEQ => 11,
    MERGE_LEN => 12,
    FAIDX_PADDING => 13,
    GEN_REF_1 => 14,
    GEN_REF_2 => 15,
    #-------------
    MICROHOM_LEN => 16,
    MICROHOM_MATCH => 17,
    JXN_BASES => 18,
    SV_SIZE => 19 ,
    #-------------
    N_TOTAL => 20,
    N_SPLITS => 21,
    N_GAPS => 22,
    N_OUTER_CLIPS => 23,
    N_DUPLEX => 24,
    N_DUPLEX_ALL => 25,
    NET_STRAND_COUNT => 26,
    NET_STRAND_COUNT_ALL => 27,
    N_SHARED_PROPER => 28,
    N_SHARED_PROPER_ALL => 29,
    #------------- 
    IS_MERGED => 30,
    SEQ_LEN => 31,
    UMI => 32,
    MAPQ => 33,
    #-------------
    CHUNK_OFFSET => 34,
    CHUNK_SIZE => 35,
    #-------------
    SAMPLE => 36
};

# working variables
my ($prevChromSide, $prevPos1, $groupIndex) = ("", 0, 0);

# loop pre-sorted SVs and commit as chrom-strand-specific proximity groups
# the code forces a break between SVs that could not possibly match at POS_1
while(my $line = <STDIN>){
    chomp $line;
    my @f = split("\t", $line);
    my $chromSide = join(":", @f[CHROM_2, SIDE_2, CHROM_1, SIDE_1]);
    
    # break conditions
    if(!$prevChromSide or
        $prevChromSide ne $chromSide or
        $f[POS_1] - $prevPos1 > $MAX_MAX_TLEN){
        $groupIndex++;
    }
    
    # print line and get ready for the next
    print join("\t", 
        $line, 
        0,   # N_MATCHES
        ".", # MATCHING_SVS
        "$f[SAMPLE]:$f[SV_ID]", # svKey
        $groupIndex # groupIndex
    ), "\n";    
    $prevChromSide = $chromSide;
    $prevPos1 = $f[POS_1]
}
