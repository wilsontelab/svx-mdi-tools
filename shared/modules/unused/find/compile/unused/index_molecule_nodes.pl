use strict;
use warnings;

# created a sorted nodes file with an index to rapidly retrieve
# all nodes claimed by a given molecule

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
my ($prevMolId, $offset, $molChunkSize) = ("", 0, 0);

# output file handles
my $outFile = "$ENV{COMPILE_PREFIX}.nodes_by_molecule.txt";
my $idxFile = join(".", $outFile, 'index');
open my $outH, ">", $outFile or die "could not open $outFile for writing\n";
open my $idxH, ">", $idxFile or die "could not open $idxFile for writing\n";
print $idxH join("\t", qw(molecule offset size)), "\n";

# loop pre-sorted nodes and commit as molecule groups
while(my $line = <STDIN>){
    print $outH $line;    
    my @f = split("\t", $line, MOL_ID + 2);
    if($prevMolId and $prevMolId != $f[MOL_ID]){
        print $idxH join("\t", $prevMolId, $offset, $molChunkSize), "\n";
        $offset += $molChunkSize;        
        $molChunkSize = 0;
    }
    $prevMolId = $f[MOL_ID];
    $molChunkSize += length($line);
}

# commit the last molecule
print $idxH join("\t", $prevMolId, $offset, $molChunkSize), "\n";

# finish up
close $outH;
close $idxH;
