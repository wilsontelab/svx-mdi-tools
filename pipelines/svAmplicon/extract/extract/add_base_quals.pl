use strict;
use warnings;

# calculate the average base quality for each alignment segment

# initialize reporting
our $script = "add_base_quals";
our $error  = "$script error";
my ($nMolecules) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# constants
use constant {
    MOL_ID => 0, # amplicon fields
    CHROM1 => 1,
    STRAND1 => 2,
    POS1 => 3,
    SEQ1 => 4,
    QUAL1 => 5,
    READ_N1 => 6,
    CHROM2 => 7,
    STRAND2 => 8,
    POS2 => 9,
    SEQ2 => 10,
    QUAL2 => 11,
    READ_N2 => 12,
    MAPQ => 13,
    CIGAR => 14,
    ALNQ => 15,
    EDGE_TYPE => 16,
    EVENT_SIZE => 17,
    INSERT_SIZE => 18,
    JXN_BASES => 19,
    AMPLICON_ID => 20,
    MERGE_LEVEL => 21,
    OVERLAP => 22,
    IS_REFERENCE => 23,
    N_READ_PAIRS => 24,
    BASE_QUAL => 25, # added by us to alignment edges only
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    MERGE_FAILURE => "M",
    PROPER        => "P",
    REJECTED_INDEL => "R"
};

# process data by molecule over multiple parallel threads
my ($prevQBases, $prevMolId, @lines) = (0);
my $lineN = 0;
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line);

    # commit molecules
    my $molId = $line[AMPLICON_ID].":".$line[MOL_ID];
    if($prevMolId and $prevMolId ne $molId){
        $nMolecules++;
        print join("", map { join("\t", @$_)."\n" } @lines);
        @lines = ();
        $prevQBases = 0;
    } 

    # assign alignment base qualities
    if($line[EDGE_TYPE] eq ALIGNMENT){ # READ_N1/2, QUAL1/2 are identical between alignment nodes
        if(@lines >= 2 and $line[READ_N2] != $lines[$#lines - 1][READ_N2]){ # must be READ_N2, in case previous alignment was a merge failure
            $prevQBases = 0; # reset qryPos on read2 start
        }
        if($prevQBases == 0){ # handle leading S on first aln only; splits do NOT have leading S operation, CIGAR SVs do
            if($line[STRAND2] eq "+"){
                $line[CIGAR] =~ m/^(\d+)S/ and $prevQBases += $1;
            } else {
                $line[CIGAR] =~ m/(\d+)S$/ and $prevQBases += $1;
            }
        }
        my $nQryBases = 0;
        while ($line[CIGAR] =~ (m/(\d+)(\w)/g)) {
            ($2 eq "M" or $2 eq "I") and $nQryBases += $1;
        }
        my $qual = substr($line[QUAL2], $prevQBases, $nQryBases);
        $line[BASE_QUAL] = int(getAvgQual($qual) + 0.5); # QUAL1 already reversed on "-" strand
        $prevQBases += $nQryBases; # accounts for leading bases for both split and CIGAR SVs (see above)
    } else {
        $line[BASE_QUAL] = "NA";
        $prevQBases += ($line[INSERT_SIZE] eq "NA" ? 0 : $line[INSERT_SIZE]); # merge failures have NA for size and insSize; will reset on next alignment
    }
    $prevMolId = $molId;
    push @lines, \@line;
}
$nMolecules++;
print join("", map { join("\t", @$_)."\n" } @lines);

# print summary information
printCount($nMolecules, 'nMolecules', 'input molecules');

1;
