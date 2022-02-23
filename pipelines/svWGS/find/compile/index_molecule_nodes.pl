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
    MOL_ID => 12,  # molecule-level information
    MOL_CLASS => 13,
    MOL_STRAND => 14,
    IS_OUTER_CLIP1 => 15,
    IS_OUTER_CLIP2 => 16 
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
