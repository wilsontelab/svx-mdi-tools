use strict;
use warnings;

# group read pairs into inferred source molecules
# when requested, use fuzzy logic on positions to allow for indels in outer clipped bases

# initialize reporting
our $script = "group_reads";
my ($nInputReadPairs, $nOutputGroups, $nOutputReadPairs) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
resetCountFile();

# environment variables
fillEnvVar(\my $ADJACENCY_LIMIT,  'ADJACENCY_LIMIT');

# constants
use constant {
    MOL_KEY => 0, # input fields, as output by parse_BWA.pl
    GROUP_POS1 => 1,
    GROUP_POS2 => 2, 
    SEQ1 => 3,
    QUAL1 => 4,    
    SEQ2 => 5,
    QUAL2 => 6,
    QNAME => 7,
    MOL_STRAND => 8,
    MERGE_LEVEL => 9,
    RANDOM_SORT => 10,
    #-----------
    STRAND1 => 0,
    STRAND2 => 1,
    #-----------
    UMI1 => 0, # indices within MOL_KEY
    RNAME1 => 1,
    SIDE1 => 2,
    IS_SV_CLIP1 => 3,
    UMI2 => 4, 
    RNAME2 => 5,
    SIDE2 => 6,
    IS_SV_CLIP2 => 7,
    GROUPING_STRAND => 8
};

# working variables
my ($prevKey, $prevPos1, @pos1Group);

# run the read pairs
# NB: cannot multi-thread when using multiple output lines for a molecule group, which must stay together
# not worth changing as this is not a particularly slow step in the overall pipeline
while(my $line = <STDIN>){ # expects input sorted by MOL_KEY+GROUP_POS1 (will sort for GROUP_POS2 etc. below)
    $nInputReadPairs++;
    chomp $line;
    my @pos1Pair = split("\t", $line);
    
    # split the molecule key groups by GROUP_POS1
    if($prevKey and
       ($prevKey ne $pos1Pair[MOL_KEY] or
        $prevPos1 < $pos1Pair[GROUP_POS1] - $ADJACENCY_LIMIT)){
        processPos1Group();
        @pos1Group = (); # prepare for next MOL_KEY+GROUP_POS1 group
    }
    
    # add read-pair to MOL_KEY+GROUP_POS1 fuzzy logic set
    push @pos1Group, \@pos1Pair;
    $prevKey  = $pos1Pair[MOL_KEY];    
    $prevPos1 = $pos1Pair[GROUP_POS1];
}
processPos1Group(); # finish last data

# print summary information
printCount($nInputReadPairs,  'nInputReadPairs',  'input read pairs (one per line)');
printCount($nOutputGroups,    'nOutputGroups',    'output molecule groups');
printCount($nOutputReadPairs, 'nOutputReadPairs', 'read pairs in the output molecule groups');

# process a group of readPairs with same umis,chroms,strands and nearby GROUP_POS1
sub processPos1Group {
    
    # sort the pos1 group by pos2
    # makes 123/546 and 124/546 adjacent in new list, despite presence of 123/987
    @pos1Group > 1 and @pos1Group = sort {$$a[GROUP_POS2] <=> $$b[GROUP_POS2]} @pos1Group;

    # split the GROUP_POS1 group by GROUP_POS2
    my ($prevPos2, @pos2Group, @strandCount);
    foreach my $pos2Pair(@pos1Group){
        if($prevPos2 and $prevPos2 < $$pos2Pair[GROUP_POS2] - $ADJACENCY_LIMIT){
            processPos2Group($prevPos2, \@pos2Group, \@strandCount);
            @pos2Group = (); # prepare for next GROUP_POS2 group
            @strandCount = ();
        }
        
        # add read-pair to the GROUP_POS2 fuzzy logic set
        push @pos2Group, $pos2Pair;
        $prevPos2 = $$pos2Pair[GROUP_POS2];
        $strandCount[$$pos2Pair[MOL_STRAND]]++; 
    }
    processPos2Group($prevPos2, \@pos2Group, \@strandCount); # finish last data
}

# process a group of readPairs with the same umis,chroms,strands and nearby GROUP_POS1+GROUP_POS2
# i.e., an inferred source DNA molecule
sub processPos2Group {
    my ($prevPos2, $pos2Group, $strandCount) = @_;

    # sort the pos2 group by strand of the source molecule
    @$pos2Group > 1 and @$pos2Group = sort {
        $$a[MOL_STRAND]  <=> $$b[MOL_STRAND] or
        $$b[MERGE_LEVEL] <=> $$a[MERGE_LEVEL] or # merged molecules come first on a strand
        $$a[RANDOM_SORT] <=> $$b[RANDOM_SORT] # but otherwise are randomized for possible downsampling
    } @$pos2Group;     

    # get molecule level information
    my @molKey = split(",", $$pos2Group[0][MOL_KEY], SIDE2);
  
    # assign a unique identifier to the source molecule
    $nOutputGroups++;

    # print output, one read pair per line, sorted by molecule+strand
    foreach my $readPair(@$pos2Group){
        $nOutputReadPairs++;
        print join("\t", $nOutputReadPairs,
                         @$readPair[MOL_STRAND,
                                    SEQ1, QUAL1,
                                    SEQ2, QUAL2,
                                    QNAME]), "\n";
    }
    
    # print a trailing line with group-level information required for consensus making
    print join("\t", '@M',
        $nOutputGroups,
        @molKey[UMI1, UMI2], # used after re-mapping during endpoint uniqueness assessment (could be 1,1 if no UMIs)
        $$strandCount[STRAND1] || 0,
        $$strandCount[STRAND2] || 0
    ), "\n";
}
