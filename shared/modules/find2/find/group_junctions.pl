use strict;
use warnings;

# break junctions into smaller continuity groups prior to aggregating junction evidence
# continue to apply SV filters as possible at this stage

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);

# environment variables
fillEnvVar(\my $MAX_TLENS,    'MAX_TLENS');
fillEnvVar(\my $MIN_COVERAGE, 'MIN_COVERAGE');

# parse the group splitting threshold as the largest library fragment size
my @MAX_TLENS = sort { $b <=> $a } split(/\s/, $MAX_TLENS);
my $MAX_MAX_TLEN = $MAX_TLENS[0];

# working variables
my ($prevChromSide, $prevPos1, $groupIndex, @p1Group) = ("", 0, 0);

# loop pre-sorted junctions and commit as chrom-strand-specific proximity groups
# this first loop forces a break between junctions that could not possibly match at POS_1
while(my $longLine = <STDIN>){
    chomp $longLine;
    my ($chromSide, $pos1, $pos2, $shortLine) = split("\t", $longLine, 4);
    if($prevChromSide and
        ($prevChromSide ne $chromSide or
         $pos1 - $prevPos1 > $MAX_MAX_TLEN)){
        processP1Group();
        @p1Group = ();
    }
    push @p1Group, [$pos2, $shortLine];
    $prevChromSide = $chromSide;
    $prevPos1 = $pos1;
}
processP1Group(); # finish the last group

# this second loop forces a break between junctions that could not possibly match at POS_2
sub processP1Group {
    @p1Group == 1 and return printContinuityGroup([$p1Group[0][1]]);
    my ($prevPos2, @p2Group) = (0);
    foreach my $jxn(sort { $$a[0] <=> $$b[0] } @p1Group){
        if($prevPos2 and
            $$jxn[0] - $prevPos2 > $MAX_MAX_TLEN){
            printContinuityGroup(\@p2Group);
            @p2Group = ();
        }          
        push @p2Group, $$jxn[1];
        $prevPos2 = $$jxn[0];
    }
    printContinuityGroup(\@p2Group);
}

# print the original junction lines plus a group index for use in R
sub printContinuityGroup {
    my ($p2Group) = @_;
    @$p2Group >= $MIN_COVERAGE or return; # group cannot possibly have an SV with sufficient coverage
    $groupIndex++; # used for parallel SV finding in call_svs.R
    foreach my $shortLine(@$p2Group){
        print join("\t", $shortLine, $groupIndex), "\n";          
    }
}
