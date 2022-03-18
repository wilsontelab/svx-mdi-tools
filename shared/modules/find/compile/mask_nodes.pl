use strict;
use warnings;

# fill the SHARED_PROPER evidence column with the number
# of endpoints the SV molecule shared with proper molecules

# initialize reporting
our $script = "mask_nodes";
our $error  = "$script error";

# constants
use constant {
    NODE => 0, # node-level data
    CLIP_LEN => 1,
    CLIP_SEQ => 2,
    NODE_CLASS => 3,
    PARTNER => 4, # chromI:side for the other node in a junction or molecule
    #---------------
    JXN_TYPE => 5, # edge/junction-level data
    JXN_N => 6,
    #---------------
    FLAG => 7, # alignment-level data
    POS => 8,
    MAPQ => 9,
    CIGAR => 10,
    SEQ => 11,
    ALN_N => 12,
    #---------------
    MOL_ID => 13, # molecule-level information  
    UMI => 14,
    IS_MERGED => 15,
    IS_DUPLEX => 16,
    STRAND_COUNT1 => 17,
    STRAND_COUNT2 => 18,
    MOL_CLASS => 19,
    MOL_STRAND => 20,
    IS_OUTER_CLIP1 => 21,
    IS_OUTER_CLIP2 => 22,
    TARGET_CLASS => 23,
    SHARED_PROPER => 24,
    #---------------
    OUT_POS_1 => 25,
    OUT_POS_2 => 26,
};

# working variables
our (%knownEndpoints);

# load the list of junctions with ends that match proper pairs
my $mpFile = "$ENV{COMPILE_PREFIX}.matchedProper.gz";
open my $endpointH, "-|", "gunzip -c $mpFile" or die "$error: could not open $mpFile: $!\n";
while(my $line = <$endpointH>){
    chomp $line;
    my ($molId, $nSharedProper) = split("\t", $line);
    $knownEndpoints{$molId} = $nSharedProper;
}
close $endpointH;

# fill the SHARED_PROPER column of a node stream
while(my $line = <STDIN>){
    my @node = split("\t", $line);
    $node[SHARED_PROPER] = $knownEndpoints{$node[MOL_ID]} || 0;
    print join("\t", @node);
}
