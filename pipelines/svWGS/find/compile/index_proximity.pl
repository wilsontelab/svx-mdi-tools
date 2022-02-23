use strict;
use warnings;

# created a sorted nodes file with an index to rapidly retrieve
# all nodes within a coordinate proximity of query node
# the returned nodes will:
#   be in the same orientation as the query node
#   include nodes to the left and right of the query node
#   include the query node itself

# constants
use constant {
    CHROM_STRAND => 0,
    POSITION => 1,
    #---------------
    NODE => 2, # node-level data
    CLIP_LEN => 3,
    CLIP_SEQ => 4,
    NODE_CLASS => 5,
    #---------------
    JXN_TYPE => 6, # edge/junction-level data
    JXN_N => 7,
    #---------------
    FLAG => 8, # alignment-level data
    POS => 9,
    MAPQ => 10,
    CIGAR => 11,
    SEQ => 12,
    ALN_N => 13,
    #---------------
    MOL_ID => 14,  # molecule-level information
    MOL_CLASS => 15,
    MOL_STRAND => 16,
    IS_OUTER_CLIP1 => 17,
    IS_OUTER_CLIP2 => 18,
    #---------------
    INITIALIZED => 0,
    PENDING => 1,
    COMMITTED => 2
};

# operating parameters
my $maxTLen = $ENV{MAX_TLEN};
#my $minAlnLen = $ENV{MIN_ALN_LEN} || 30; # estimate of the smallest number of bases BWA needs to align a base segment
my $separationLimit = $maxTLen * 1.1; # - 2 * $minAlnLen;

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
    my $outLine = join("\t", @f[NODE..IS_OUTER_CLIP2])."\n";
    print $outH $outLine;
    
    # break condition between nodes
    if($prevChromStrand and ( 
            $prevChromStrand ne $f[CHROM_STRAND] or
            abs($f[POSITION] - $nodes[0]{position}) > $separationLimit
        )){
        processProximityGroups($f[POSITION]); # do work on the current node set
        if($prevChromStrand ne $f[CHROM_STRAND]){ # start a new chrom+side afresh
            @nodes = ();
            %committed = ();
        } else {
            while(@nodes > 0 and $nodes[0]{status} == COMMITTED){ # remove nodes we have fully characterized
                shift @nodes;                                     # keeps memory small
            }            
        }
    }
    
    # add the current node to the growing set
    $prevChromStrand = $f[CHROM_STRAND];
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
