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
our $error  = "$script error";
my ($nJxnNodes, $nRepeatedNodes) = (0) x 20;

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/../_workflow/workflow.pl";
require "$scriptDir/../common/utilities.pl";
resetCountFile();

# constants
use constant {
    SV_ID => 0, # SV-level data
    NODE_N => 1, 
    IS_JUNCTION_NODE => 2, 
    IS_SEED_NODE => 3,
    IS_REF_NODE => 4, 
    IS_RC => 5, 
    IS_REPEAT => 6, 
    N_COLLAPSED => 7,
    #---------------   
    NODE => 8, # node-level data
    CLIP_LEN => 9,
    CLIP_SEQ => 10,
    NODE_CLASS => 11,
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
    MOL_ID => 20,  # molecule-level information
    MOL_CLASS => 21,
    MOL_STRAND => 22,
    IS_OUTER_CLIP1 => 23,
    IS_OUTER_CLIP2 => 24 
};

# working variables
my (%seen);
my $lastSplitI = MOL_ID + 2;

# thread through nodes, which come to us in SV_ID order
while (my $line = <STDIN>) {
    my @node = split("\t", $line, $lastSplitI); # retains newline
    if ($node[IS_JUNCTION_NODE]){ # NB: only junction nodes are marked as repeat/no-repeat
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
    } else {
        print $line;
    }
}

printCount($nJxnNodes,      'nJxnNodes',      'input junction node entries');
printCount($nRepeatedNodes, 'nRepeatedNodes', 'repeated junction node entries');
