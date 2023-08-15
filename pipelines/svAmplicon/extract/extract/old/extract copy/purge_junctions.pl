use strict;
use warnings;

# purge merge-failure and too-small indels from the node stream
# by collapsing the alignments that flank these unusuable/unwanted events

# along the way, calculate the average base quality for each alignment segment

# initialize reporting
our $script = "purge_junctions";
our $error  = "$script error";
my ($nMolecules) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $MIN_SV_SIZE, 'MIN_SV_SIZE');

# constants
use constant {
    MOL_ID => 0, # amplicon fields
    CHROM1 => 1,
    STRAND1 => 2,
    POS1 => 3,
    SEQ1 => 4,
    QUAL1 => 5,
    CIGAR1 => 6,
    READ_N1 => 7,
    CHROM2 => 8,
    STRAND2 => 9,
    POS2 => 10,
    SEQ2 => 11,
    QUAL2 => 12,
    CIGAR2 => 13,
    READ_N2 => 14,
    EDGE_TYPE => 15,
    MAPQ => 16,
    SIZE => 17,
    INS_SIZE => 18,
    AMPLICON_ID => 19,
    MERGE_LEVEL => 20,
    OVERLAP => 21,
    IS_REF => 22,
    N_READ_PAIRS => 23,
    BASE_QUAL => 24, # added by us
    # TODO: could add ALN_SEQ, following BASE_QUAL calc below, if it would be useful to pass to app
    # however, final moleculeTypes table will have the canonical strand read1 (and maybe read2) present
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
    FUSED_MERGE_FAILURE => "F",
    REJECTED_INDEL => "R",
    FUSED_MERGE_FAILURE_REJECTED_INDEL => "Q"
};

# process data by molecule over multiple parallel threads
my ($prevQBases, $prevMolId, @lines) = (0);
while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line);
    # print join("\t", @line[MOL_ID..POS1], @line[CIGAR1..POS2], @line[CIGAR2..N_READ_PAIRS]), "\n"; # FOR DEBUGGING
    # next;

    # commit molecules
    my $molId = $line[AMPLICON_ID].":".$line[MOL_ID];
    if($prevMolId and $prevMolId ne $molId){
        $nMolecules++;
        print join("", map { join("\t", @$_)."\n" } @lines);
        @lines = ();
        $prevQBases = 0;
    } 

    # assign alignment base qualities
    if($line[EDGE_TYPE] eq ALIGNMENT){ # CIGAR1/2, READ_N1/2, QUAL1/2 are identical between alignment nodes
        if(@lines >= 2 and $line[READ_N2] != $lines[$#lines - 1][READ_N2]){ # must be READ_N2, in case previous alignment was a merge failure
            $prevQBases = 0; # reset qryPos on read2
        }
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
        my $qual = substr($line[QUAL2], $prevQBases, $nQryBases);
        $line[BASE_QUAL] = int(getAvgQual($qual) + 0.5); # QUAL1 already reversed on "-" strand
        $prevQBases += $nQryBases; # accounts for leading bases for both split and CIGAR SVs (see above)
    } else {
        $line[BASE_QUAL] = "NA";
        $prevQBases += ($line[INS_SIZE] eq "NA" ? 0 : $line[INS_SIZE]); # merge failures have NA for size and insSize; will resent on next alignment
    }
    # print join("\t", @line[MOL_ID..POS1], @line[CIGAR1..POS2], @line[CIGAR2..BASE_QUAL]), "\n"; # FOR DEBUGGING
    # next;

    # collapse segments
    if(@lines >= 2){
        my $prevLine = $lines[$#lines];
        my $prevPrevLine = $lines[$#lines - 1];

        # merge failures
        if($$prevLine[EDGE_TYPE] eq MERGE_FAILURE){
            $line[CHROM1]   = $$prevPrevLine[CHROM1];
            $line[STRAND1]  = $$prevPrevLine[STRAND1];
            $line[POS1]     = $$prevPrevLine[POS1];
            $line[SEQ1]     = $$prevPrevLine[SEQ1];
            $line[QUAL1]    = $$prevPrevLine[QUAL1];
            $line[CIGAR1]   = $$prevPrevLine[CIGAR1];
            $line[READ_N1]  = $$prevPrevLine[READ_N1];
            $line[EDGE_TYPE] = $$prevPrevLine[EDGE_TYPE] eq REJECTED_INDEL ? FUSED_MERGE_FAILURE_REJECTED_INDEL : FUSED_MERGE_FAILURE;
            $line[MAPQ] = $$prevLine[MAPQ];
            $line[SIZE] = $$prevPrevLine[SIZE] + $line[SIZE] - ($line[OVERLAP] eq "NA" ? 0 : $line[OVERLAP]);
            $line[INS_SIZE] = 0;
            $line[BASE_QUAL] = min($$prevPrevLine[BASE_QUAL], $line[BASE_QUAL]);
            @lines = @lines == 2 ? () : @lines[0..$#lines - 2];

        # ~matched size ins+del
        } elsif(($$prevLine[EDGE_TYPE] eq DELETION or $$prevLine[EDGE_TYPE] eq INSERTION) and 
                 abs($$prevLine[SIZE] - $$prevLine[INS_SIZE]) <= $MIN_SV_SIZE){
            $$prevPrevLine[CIGAR2] =~ s/\d+S$//g;
            $line[CIGAR1] =~ s/^\d+S//g;
            my $D = $$prevLine[INS_SIZE] < 0 ? $$prevLine[SIZE] + $$prevLine[INS_SIZE] : $$prevLine[SIZE];
            my $I = $$prevLine[INS_SIZE] < 0 ? 0 : $$prevLine[INS_SIZE];
            my $cigar = join("", 
                $$prevPrevLine[CIGAR2], 
                ($D ? $D."D": ()), 
                ($I ? $I."I": ()), 
                $line[CIGAR1]
            );
            $line[POS1]   = $$prevPrevLine[POS1];
            $line[CIGAR1] = $$prevPrevLine[EDGE_TYPE] eq FUSED_MERGE_FAILURE ? $$prevPrevLine[CIGAR1] : $cigar;
            $line[CIGAR2] = $cigar;
            $line[EDGE_TYPE] = $$prevPrevLine[EDGE_TYPE] eq FUSED_MERGE_FAILURE ? FUSED_MERGE_FAILURE_REJECTED_INDEL : REJECTED_INDEL;
            $line[MAPQ]   = $$prevLine[MAPQ];
            $line[SIZE]   = $$prevPrevLine[SIZE] + $$prevLine[SIZE] + $line[SIZE];         
            $line[INS_SIZE] = 0;
            $line[BASE_QUAL] = min($$prevPrevLine[BASE_QUAL], $line[BASE_QUAL]);
            @lines = @lines == 2 ? () : @lines[0..$#lines - 2];
        } 
    }
    $prevMolId = $molId;
    push @lines, \@line;
}
$nMolecules++;
print join("", map { join("\t", @$_)."\n" } @lines);

# print summary information
printCount($nMolecules, 'nMolecules', 'input molecules');

1;
