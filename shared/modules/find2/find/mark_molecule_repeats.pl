use strict;
use warnings;

# add a mark to SV nodes that were also claimed by an earlier numbered SV

# thus, IS_REPEAT == 0 yields a table where junction instances are only ever present once
# two ways this repeating junction pattern can happen are:
#   1) a very long molecule fails proximity match to a very short one (but each match intermediate length molecules)
#   2) an alignment of one split molecule instance is aberrant enough to escape PURGE_DISTANCE aggregation
# junction instances will first get added to the "best" query junctions (most frequent, splits)
# also, a junction can never be both IS_REPEAT and IS_SEED_JUNCTION; thus, later redundant SVs will have different seeds

# initialize reporting
our $script = "mark_molecule_repeats";
my ($nJxnNodes, $nRepeatedNodes) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
resetCountFile();

# constants
use constant {
    SV_ID => 0, # SV-level data
    NODE_N => 1, 
    IS_SEED_NODE => 2,
    IS_REF_NODE => 3, 
    IS_RC => 4, 
    IS_REPEAT => 5, 
    N_COLLAPSED => 6,
    #---------------
    NODE => 7, # node-level data
    CLIP_LEN => 8,
    CLIP_SEQ => 9,
    NODE_CLASS => 10,
    PARTNER => 11, # chromI:side for the other node in a junction or molecule
    #---------------
    JXN_TYPE => 12, # edge/junction-level data
    JXN_N => 13,
    #---------------
    FLAG => 14, # alignment-level data
    POS => 15,
    MAPQ => 16,
    CIGAR => 17,
    SEQ => 18,
    ALN_N => 19,
    #---------------
    MOL_ID => 20, # molecule-level information  
    UMI => 21,
    IS_MERGED => 22,
    IS_DUPLEX => 23,
    STRAND_COUNT1 => 24,
    STRAND_COUNT2 => 25,
    MOL_CLASS => 26,
    MOL_STRAND => 27,
    IS_OUTER_CLIP1 => 28,
    IS_OUTER_CLIP2 => 29,
    TARGET_CLASS => 30,
    SHARED_PROPER => 31,
    #---------------
    OUT_POS_1 => 32,
    OUT_POS_2 => 33
};

# working variables
my (%seen);
my $lastSplitI = MOL_ID + 2;

# thread through nodes, which come to us in SV_ID order
while (my $line = <STDIN>) {
    my @node = split("\t", $line, $lastSplitI); # retains newline
    $nJxnNodes++;
    my $alnId = join(":", @node[MOL_ID, JXN_N, ALN_N]);
    if ($seen{$alnId}) {
        $nRepeatedNodes++;
        $node[IS_REPEAT] = 1; # NB: does/must not change the character length of the line
        print join("\t", @node);
    } else {
        $seen{$alnId} = 1;
        print $line;
    }
}

printCount($nJxnNodes,      'nJxnNodes',      'input junction node entries');
printCount($nRepeatedNodes, 'nRepeatedNodes', 'repeated junction node entries');
