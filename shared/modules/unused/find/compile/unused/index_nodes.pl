use strict;
use warnings;

# created a sorted nodes file with an index to rapidly retrieve
# all molecules that claim a given node

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
    UMI => 12, # molecule-level information
    MOL_ID => 13,       
    IS_DUPLEX => 14,     
    STRAND_COUNT1 => 15,
    STRAND_COUNT2 => 16,
    IS_MERGED => 17,
    TARGET_CLASS => 18,
    MOL_CLASS => 19,
    SHARED_PROPER => 20,
    IS_OUTER_CLIP1 => 21,
    IS_OUTER_CLIP2 => 22,
    IS_ORPHAN => 23,
    MOL_STRAND => 24
};

# working variables
my ($prevNodeName, $offset, $nodeChunkSize) = ("", 0, 0);

# output file handles
my $outFile = "$ENV{COMPILE_PREFIX}.nodes_by_node.txt";
my $indexFile = join(".", $outFile, 'index');
open my $outH, ">", $outFile or die "could not open $outFile for writing\n";
open my $idxH, ">", $indexFile or die "could not open $indexFile for writing\n";
print $idxH join("\t", qw(node offset size)), "\n";

# loop pre-sorted nodes and commit as node groups
while(my $line = <STDIN>){
    print $outH $line;    
    my @f = split("\t", $line, NODE + 2);
    if($prevNodeName and $prevNodeName ne $f[NODE]){
        print $idxH join("\t", $prevNodeName, $offset, $nodeChunkSize), "\n";
        $offset += $nodeChunkSize;        
        $nodeChunkSize = 0;
    }
    $prevNodeName = $f[NODE];
    $nodeChunkSize += length($line);
}

# commit the last node
print $idxH join("\t", $prevNodeName, $offset, $nodeChunkSize), "\n";

# finish up
close $outH;
close $idxH;

