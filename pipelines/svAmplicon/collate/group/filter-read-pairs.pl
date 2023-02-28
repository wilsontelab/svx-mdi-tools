use strict;
use warnings;

# apply fuzzy positions to keep read pairs matching an amplicon index
# finish the process of repairing the outer ends to match the amplicon
# (outer ends are never informative for a PCR product, so no information is lost here)

# environment variables
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
fillEnvVar(\our $COLLATE_PREFIX,   'COLLATE_PREFIX');
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svAmplicon";
map { require "$perlUtilDir/$_.pl" } qw(amplicons);
use vars qw(@amplicons);

# constants
use constant {
    CHROM1  => 0, # shared leader columns
    SIDE1   => 1,
    POS1    => 2,
    CHROM2  => 3,
    SIDE2   => 4,
    POS2    => 5,
    #--------------
    COUNT   => 6, # ALLOWED_FILE extended columns
    AMPLICON=> 7,
    INDEX   => 8,
    PATCH1  => 9,
    PATCH2  => 10,
    #--------------
    UMI1    => 6, #INTERIM_FILE extended columns
    UMI2    => 7,
    SEQ1    => 8,
    QUAL1   => 9,
    SEQ2    => 10,
    QUAL2   => 11,
    OVERLAP => 12,
    MERGED  => 13,
    #-------------
    AMP_AMPLICON_ID => 0, # amplicon fields
    AMP_PROPER => 1,
    AMP_MOL_COUNT => 2,
    AMP_CHROM1 => 3,
    AMP_SIDE1 => 4,
    AMP_POS1 => 5,
    AMP_REF1 => 6,
    AMP_CHROM2 => 7,
    AMP_SIDE2 => 8,
    AMP_POS2 => 9,
    AMP_REF2 => 10,
};

# collect the list of allowed outer node pairs
my %allowedNodePairs;
open my $inH, "<", $ENV{ALLOWED_FILE} or die "not found: $ENV{ALLOWED_FILE}\n";
while(my $nodePair = <$inH>){
    chomp $nodePair;
    my @nodePair = split("\t", $nodePair);
    my $key = join("\t", @nodePair[CHROM1..POS2]);
    $allowedNodePairs{$key} = \@nodePair;
}
close $inH;

# filter the INTERIM_FILE molecules stream
sub applyEndPatch {
    my ($patch, $mol, $side, $seq, $qual) = @_;
    $patch eq "0" and return;
    if(substr($patch, 0, 1) eq "-"){ # bases to trim
        my $trim = $patch + 0; # a negative integer
        if($$mol[$side] eq "R"){
            $$mol[$seq]  = substr($$mol[$seq],  -$trim);
            $$mol[$qual] = substr($$mol[$qual], -$trim);
        } else {
            $$mol[$seq]  = substr($$mol[$seq],  0, $trim);
            $$mol[$qual] = substr($$mol[$qual], 0, $trim);
        }
    } else { # padding bases to add
        my $qualPatch = "#" x length($patch);
        if($$mol[$side] eq "R"){
            $$mol[$seq]  = "$patch$$mol[$seq]";
            $$mol[$qual] = "$qualPatch$$mol[$qual]";
        } else {
            $$mol[$seq]  = "$$mol[$seq]$patch";
            $$mol[$qual] = "$$mol[$qual]$qualPatch";
        }
    }
}
while(my $mol = <STDIN>){
    chomp $mol;
    my @mol = split("\t", $mol);
    my $key = join("\t", @mol[CHROM1..POS2]);
    my $nodePair = $allowedNodePairs{$key};
    !$nodePair and next; # apply the molecule filter
    applyEndPatch($$nodePair[PATCH1], \@mol, SIDE1, SEQ1, QUAL1); # patch the molecule ends
    if($mol[SEQ2] eq "*"){
        applyEndPatch($$nodePair[PATCH2], \@mol, SIDE2, SEQ1, QUAL1);
    } else {
        applyEndPatch($$nodePair[PATCH2], \@mol, SIDE2, SEQ2, QUAL2);  
    }
    my $amplicon = $amplicons[$$nodePair[AMPLICON]];
    my $isReference = 0; # notPossible stays here, of course, no read pair can be reference
    if($$amplicon[AMP_PROPER] eq "expectOverlap"){ # unmerged expectOverlap reads cannot possibly have the reference sequence (they would have merged)
        $mol[SEQ1] eq $$amplicon[AMP_REF1] and $mol[SEQ2] eq "*" and $isReference = 1;
    } elsif($$amplicon[AMP_PROPER] eq "expectGaps") {
        $mol[SEQ1] eq $$amplicon[AMP_REF1] and $mol[SEQ2] eq $$amplicon[AMP_REF2] and $isReference = 1;
    }
    print join("\t", 
        $$nodePair[AMPLICON],     
        @mol[SEQ1, SEQ2, QUAL1, QUAL2, OVERLAP],
        $isReference
    )."\n";
}
