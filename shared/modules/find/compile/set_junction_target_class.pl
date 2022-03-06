use strict;
use warnings;

# add a TARGET_CLASS field that applies to the two nodes that form an SV junction
# NB: this is distinct from the similarly named field that applies to the outer ends of _molecules_

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms targets);

# environment variables
fillEnvVar(\our $TARGETS_BED,      'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,   'REGION_PADDING');
fillEnvVar(\our $TARGET_SCALAR,    'TARGET_SCALAR', 1, 10); # use 10 bp target resolution for svCapture targets

# constants
use constant {
    JXN_N => 0,
    MOL_ID => 1,
    NODE1 => 2,
    NODE2 => 3,
    JXN_TYPE => 4
};

# load the target regions
loadTargetRegions('quiet');
use vars qw(%revChromIndex $nRegions);

# loop pre-sorted nodes and commit as molecule groups
while(my $line = <STDIN>){
    chomp $line;
    my @f = split("\t", $line);
    my ($chromI1, $side1, $pos1) = split(':', $f[NODE1]);
    my ($chromI2, $side2, $pos2) = split(':', $f[NODE2]);
    print join("\t", 
        $f[NODE1], 
        $f[NODE2], 
        $f[JXN_TYPE], 
        $nRegions ? getTargetClass($revChromIndex{$chromI1}, $revChromIndex{$chromI2}, $pos1, $pos2) : "*", 
        $f[MOL_ID]
    ), "\n";
}
