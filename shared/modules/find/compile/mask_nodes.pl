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
    #---------------
    JXN_TYPE => 4, # edge/junction-level data
    JXN_N => 5,
    #---------------
    FLAG => 6, # alignment-level data
    POS => 7,
    MAPQ => 8,
    CIGAR => 9,
    SEQ => 10,
    ALN_N => 11,
    #---------------
    MOL_ID => 12, # molecule-level information  
    UMI => 13,
    IS_MERGED => 14,
    IS_DUPLEX => 15,
    STRAND_COUNT1 => 16,
    STRAND_COUNT2 => 17,
    MOL_CLASS => 18,
    MOL_STRAND => 19,
    IS_OUTER_CLIP1 => 20,
    IS_OUTER_CLIP2 => 21,
    TARGET_CLASS => 22,
    SHARED_PROPER => 23,
    # sometimes, not always, is additional columns
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
    chomp;
    my @node = split("\t", $line);
    $node[SHARED_PROPER] = $knownEndpoints{$node[MOL_ID]} || 0;
    print join("\t", @node), "\n";
}
