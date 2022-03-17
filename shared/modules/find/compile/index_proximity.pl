use strict;
use warnings;

# create a sorted junction nodes file with an index to rapidly retrieve
# all potentially compatible nodes within coordinate proximity of a query node
# the returned nodes for a query will:
#   be in the same orientation as the query node
#   have the same partner chrom+side as the query node
#   include nodes to the left and right of the query node
#   include the query node itself
#   only include junction nodes (i.e, splits and gaps, not outer clips)

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
    CHROM_STRAND => 25,
    POSITION => 26,
    #---------------
    INITIALIZED => 0,
    PENDING => 1,
    COMMITTED => 2
};

# operating parameters
my $maxTLen = $ENV{MAX_TLEN};
my $separationLimit = $maxTLen * 1.1;

# working variables
my ($prevChromStrand, $offset, @nodes, %committed) = ("", 0);

# output file handles
my $outFile = "$ENV{COMPILE_PREFIX}.nodes_by_proximity.txt";
my $idxFile = join(".", $outFile, 'index');
open my $outH, "|-", "slurp -s 50M -o $outFile" or die "could not open $outFile for writing\n";
open my $idxH, "|-", "slurp -s 50M -o $idxFile" or die "could not open $idxFile for writing\n";
print $idxH join("\t", qw(node offset size)), "\n";

# loop pre-sorted nodes and commit as chrom-strand-specific proximity groups
# i.e., for each node, index those nodes within a distance consistent with a single junction
while(my $line = <STDIN>){
    chomp $line;
    my @f = split("\t", $line);
    my $outLine = join("\t", @f[NODE..SHARED_PROPER])."\n";
    print $outH $outLine;
    my $chromStrand = join(":", @f[PARTNER, CHROM_STRAND]);
    
    # break condition between nodes
    if($prevChromStrand and ( 
            $prevChromStrand ne $chromStrand or
            abs($f[POSITION] - $nodes[0]{position}) > $separationLimit
        )){
        processProximityGroups($f[POSITION]); # do work on the current node set
        if($prevChromStrand ne $chromStrand){ # start a new chrom+side afresh
            @nodes = ();
            %committed = ();
        } else {
            while(@nodes > 0 and $nodes[0]{status} == COMMITTED){ # remove nodes we have fully characterized
                shift @nodes;                                     # keeps memory small
            }            
        }
    }
    
    # add the current node to the growing set
    $prevChromStrand = $chromStrand;
    my $length = length($outLine);
    push @nodes, {
        name     => $f[NODE],
        offset   => $offset,
        length   => $length,
        position => $f[POSITION],
        status   => INITIALIZED,
        pending  => undef
    };
    $offset += $length;
}

# commit the last set of nodes
processProximityGroups(1e9);

# report the set of nodes within an allowable distance of a query node
# NB: the query node is included in the indexed node set
sub processProximityGroups {
    my ($currPos) = @_;
    
    # build the group to the left of each node
    foreach my $i(0..$#nodes){
        $nodes[$i]{status} == PENDING and next;        
        my $chunkSize = 0;
        foreach my $j(0..$i){ $chunkSize += $nodes[$j]{length} }        
        $nodes[$i]{pending} = [$nodes[0]{offset}, $chunkSize];
        $nodes[$i]{status} = PENDING;        
    }
    
    # add the group the right of each node and commit when ready
    foreach my $i(0..$#nodes){
        $nodes[$i]{status} == COMMITTED and next;        
        if(abs($currPos - $nodes[$i]{position}) > $separationLimit){ # 2nd break condition
            my ($offset, $chunkSize) = @{$nodes[$i]{pending}};
            if($i < $#nodes){ # the node itself was already added to the chunk
                foreach my $j(($i + 1)..$#nodes){ $chunkSize += $nodes[$j]{length} }
            }
            $committed{$nodes[$i]{name}} or # nodes instances may exist, but get only one index line
                print $idxH join("\t", $nodes[$i]{name}, $offset, $chunkSize), "\n";
            $nodes[$i]{status} = COMMITTED;
            $committed{$nodes[$i]{name}}++;
        }    
    }
}

# finish up
close $outH;
close $idxH;
