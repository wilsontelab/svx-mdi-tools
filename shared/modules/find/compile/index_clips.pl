use strict;
use warnings;

# create a sorted clips file with an index to rapidly retrieve
# all outer clips with a requested node signature, i.e., name

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
};

# working variables
my ($prevName, $offset, $chunkSize) = ("", 0, 0);

# output file handles
my $outFile = "$ENV{COMPILE_PREFIX}.outer_clips.txt";
my $idxFile = join(".", $outFile, 'index');
open my $outH, "|-", "slurp -s 50M -o $outFile" or die "could not open $outFile for writing\n";
open my $idxH, "|-", "slurp -s 50M -o $idxFile" or die "could not open $idxFile for writing\n";
print $idxH join("\t", qw(node offset size)), "\n";

# loop pre-sorted nodes and commit as chrom-strand-specific proximity groups
# i.e., for each node, index those nodes within a distance consistent with a single junction
while(my $line = <STDIN>){ 
    print $outH $line;    
    my @f = split("\t", $line, CLIP_SEQ);

    # break condition between node signatures
    if($prevName and $prevName ne $f[NODE]){
        print $idxH join("\t", $prevName, $offset, $chunkSize), "\n";
        $offset += $chunkSize;
        $chunkSize = 0;
    }
    
    # add the current node to the growing set
    $prevName = $f[NODE];
    $chunkSize += length($line);
}

# commit the last set of nodes
print $idxH join("\t", $prevName, $offset, $chunkSize), "\n";

# finish up
close $outH;
close $idxH;
