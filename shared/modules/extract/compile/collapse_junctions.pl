use strict;
use warnings;

# compress to one line for the two nodes corresponding to a single candidate junction
# update the TARGET_CLASS field to apply to the two nodes that form an SV _junction_
# NB: this is distinct from the value in node files that describes the outer the ends of _molecules_

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms targets);

# environment variables
fillEnvVar(\our $DATA_NAME,        'DATA_NAME');
fillEnvVar(\our $TARGETS_BED,      'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,   'REGION_PADDING');
fillEnvVar(\our $TARGET_SCALAR,    'TARGET_SCALAR', 1, 10); # use 10 bp target resolution for svCapture targets

# constants
use constant {
    NODE => 0, # node-level data
    CLIP_LEN => 1,
    CLIP_SEQ => 2,
    #---------------
    FLAG => 3, # alignment-level data
    POS => 4,
    MAPQ => 5,
    CIGAR => 6,
    SEQ => 7,
    ALN_N => 8,
    #---------------
    UMI => 9,
    #===============
    NODE_CLASS => 10,
    #---------------
    JXN_TYPE => 11, # edge/junction-level data
    JXN_N => 12,
    #---------------
    MOL_ID => 13, # molecule-level information  
    IS_MERGED => 14,
    IS_DUPLEX => 15,
    STRAND_COUNT1 => 16,
    STRAND_COUNT2 => 17,
    MOL_CLASS => 18,
    MOL_STRAND => 19,
    IS_OUTER_CLIP1 => 20,
    IS_OUTER_CLIP2 => 21,
    TARGET_CLASS => 22, # converted by this script from molecule-outer-end to junction target class
    SHARED_PROPER => 23,
    #---------------
    OUT_POS_1 => 24,
    OUT_POS_2 => 25,
    #---------------
    SAMPLE => 26 # added by this script to support sample-admixed SV finding
};

# load the target regions
setCanonicalChroms();
loadTargetRegions('quiet');
use vars qw(%revChromIndex $nRegions);

# loop pre-sorted nodes and commit as candidate SV junctions
while(my $line1 = <STDIN>){ # node 1
      my $line2 = <STDIN>;  # node 2 of the same junction, pre-sorted by parse_nodes.pl in proper order
    chomp $line1;
    my @f1 = split("\t", $line1);
    my @f2 = split("\t", $line2, JXN_TYPE);
    if($nRegions){
        my ($chromI1, $side1, $pos1) = split(':', $f1[NODE]);
        my ($chromI2, $side2, $pos2) = split(':', $f2[NODE]);    
        $f1[TARGET_CLASS] = getTargetClass($revChromIndex{$chromI1}, $revChromIndex{$chromI2}, $pos1, $pos2);
    }
    print join("\t", @f1, $DATA_NAME, @f2[NODE..UMI]), "\n";
}
