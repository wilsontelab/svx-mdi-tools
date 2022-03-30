use strict;
use warnings;

# create a sorted clips file with an index to rapidly retrieve
# all outer clips with a requested node signature, i.e., name

# constants
use constant {
    NODE => 0, # node-level data
    CLIP_LEN => 1,
    CLIP_SEQ => 2,
    # #---------------
    # FLAG => 3, # alignment-level data
    # POS => 4,
    # MAPQ => 5,
    # CIGAR => 6,
    # SEQ => 7,
    # ALN_N => 8,
    # #---------------
    # UMI => 9,
    # #---------------
    # NODE_CLASS => 10,
    # #---------------
    # JXN_TYPE => 11, # edge/junction-level data
    # JXN_N => 12,
    # #---------------
    # MOL_ID => 13, # molecule-level information  
    # IS_MERGED => 14,
    # IS_DUPLEX => 15,
    # STRAND_COUNT1 => 16,
    # STRAND_COUNT2 => 17,
    # MOL_CLASS => 18,
    # MOL_STRAND => 19,
    # IS_OUTER_CLIP1 => 20,
    # IS_OUTER_CLIP2 => 21,
    # TARGET_CLASS => 22,
    # SHARED_PROPER => 23,
    # #---------------
    # OUT_POS_1 => 24,
    # OUT_POS_2 => 25
};

# working variables
my ($prevKey, $offset, $chunkSize) = ("", 0, 0);

# output file handles
my $outFile = "$ENV{COMPILE_PREFIX}.outer_clips.txt";
my $idxFile = join(".", $outFile, 'index');
open my $outH, "|-", "slurp -s 50M -o $outFile" or die "could not open $outFile for writing\n";
open my $idxH, "|-", "slurp -s 50M -o $idxFile" or die "could not open $idxFile for writing\n";
print $idxH join("\t", qw(key offset size)), "\n";

# loop pre-sorted nodes and commit as chrom-strand-specific proximity groups
# i.e., for each node, index those nodes within a distance consistent with a single junction
while(my $line = <STDIN>){ 
    print $outH $line;    
    my @f = split("\t", $line, CLIP_SEQ);

    # break condition between node signatures
    if($prevKey and $prevKey ne $f[NODE]){
        print $idxH join("\t", $prevKey, $offset, $chunkSize), "\n";
        $offset += $chunkSize;
        $chunkSize = 0;
    }
    
    # add the current node to the growing set
    $prevKey = $f[NODE];
    $chunkSize += length($line);
}

# commit the last set of nodes
print $idxH join("\t", $prevKey, $offset, $chunkSize), "\n";

# finish up
close $outH;
close $idxH;
