use strict;
use warnings;

# calculate the average base quality for each alignment segment

# initialize reporting
our $script = "append_base_qual";
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
    CIGAR1 => 6,
    CHROM2 => 7,
    STRAND2 => 8,
    POS2 => 9,
    SEQ2 => 10,
    QUAL2 => 11,
    CIGAR2 => 12,
    EDGE_TYPE => 13,
    MAPQ => 14,
    SIZE => 15,
    INS_SIZE => 16,
    BASE_QUAL => 17, # added by us
    # SUB_SEQ => 18,
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    MERGE_FAILURE => "M",
    PROPER        => "P"
};

# process data by molecule over multiple parallel threads
my ($prevQBases, $prevMolId, @lines) = (0);
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line);
    # print join("\t", @line[MOL_ID..POS1], @line[CIGAR1..POS2], @line[CIGAR2..N_READ_PAIRS]), "\n"; # FOR DEBUGGING
    # next;

    # commit molecules
    my $molId = $line[MOL_ID];
    if($prevMolId and $prevMolId ne $molId){
        $nMolecules++;
        print join("", map { join("\t", @$_)."\n" } @lines);
        @lines = ();
        $prevQBases = 0;
    } 

    # assign alignment base qualities
    if($line[EDGE_TYPE] eq ALIGNMENT){ # CIGAR1/2, READ_N1/2, QUAL1/2 are identical between alignment nodes
        if($prevQBases == 0){ # handle leading S on first aln only; splits do NOT have leading S operation, CIGAR SVs do
            if($line[STRAND2] eq "+"){
                $line[CIGAR2] =~ m/^(\d+)S/ and $prevQBases += $1;
            } else {
                $line[CIGAR2] =~ m/(\d+)S$/ and $prevQBases += $1;
            }
        }
        my $nQryBases = 0;
        while ($line[CIGAR2] =~ (m/(\d+)(\w)/g)) {
            ($2 eq "M" or $2 eq "I") and $nQryBases += $1;
        }
        $line[BASE_QUAL] = int(getAvgQual(substr($line[QUAL2], $prevQBases, $nQryBases)) + 0.5); # QUAL already reversed on "-" strand        
        # $line[SUB_SEQ] = substr($line[SEQ2], $prevQBases, $nQryBases);        
        $prevQBases += $nQryBases; # accounts for leading bases for both split and CIGAR SVs (see above)
    } else {
        $line[BASE_QUAL] = "NA"; # TODO: fill these in
        # $line[SUB_SEQ] = "NA";
        $prevQBases += ($line[INS_SIZE] eq "NA" ? 0 : $line[INS_SIZE]);
    }
    # print join("\t", @line[MOL_ID..POS1], @line[CIGAR1..POS2], @line[CIGAR2..BASE_QUAL]), "\n"; # FOR DEBUGGING
    # next;

    $prevMolId = $molId;
    push @lines, \@line;
}
$nMolecules++;
print join("", map { join("\t", @$_)."\n" } @lines);

# print summary information
printCount($nMolecules, 'nMolecules', 'input molecules');

1;
