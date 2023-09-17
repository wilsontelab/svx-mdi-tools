use strict;
use warnings;

# adjust details of edges where indels were called as SVs in the CIGAR string of a single alignment
# known as inline junctions
#   add jxnBases from read sequence when insSize > 0
#   add JxnBases and adjust node positions and insertSizes and eventSizes if microhomology is present

# note: bwa mem right justifies microhomology bases

# remember, edges come in query read order, so for any inline junction either:
#   refPos1 refPos2
#   1       4
#   4       10
#   10      20
# or 
#   20     10 # NOT negative as they are for svPore, since not using genome node positions here
#   10     4
#   4      1
# inversions and translocations are irrelevant here as they are not inline junctions, i.e., they always yield split reads

# initialize reporting
our $script = "adjust_inline_junctions";
our $error  = "$script error";
my ($nMolecules, $nJunctions, $nInlineJunctions) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general faidx);
resetCountFile();

# initialize the genome
my $faFile = $ENV{GENOME_FASTA};
loadFaidx($faFile);
open our $faH, "<", $faFile; 

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
    OVERLAP => 22, # from read merging
    IS_REFERENCE => 23,
    N_READ_PAIRS => 24,
    IS_INLINE_JXN => 25,
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
my %jxnOps = map { $_ => 1 } (TRANSLOCATION, INVERSION, DUPLICATION, DELETION, INSERTION);

# process data by molecule over multiple parallel threads
my ($prevMolId, $pendingAdjustment, @lines);

while(my $line = <STDIN>){
    chomp $line;
    my @line = split("\t", $line);

    # commit molecules
    my $molId = $line[AMPLICON_ID].":".$line[MOL_ID];
    if($prevMolId and $prevMolId ne $molId){
        $nMolecules++;
        print join("", map { join("\t", @$_)."\n" } @lines);
        @lines = ();
        $pendingAdjustment = undef;
    } 

    # watch for inline SV calls
    if($jxnOps{$line[EDGE_TYPE]}){ $nJunctions++ }
    $line[IS_INLINE_JXN] = ($jxnOps{$line[EDGE_TYPE]} and $line[JXN_BASES] =~ m/^NA:/) ? 1 : 0;
    if($line[IS_INLINE_JXN]){
        $nInlineJunctions++;
        $pendingAdjustment = 1;

    # adjust inline junctions once rightmost alignment is encountered
    } elsif($pendingAdjustment){
        my $prevAln  = $#lines - 1;
        my $jxnEdge  = $#lines;
        my $jxnQryPos = $lines[$jxnEdge][JXN_BASES];
        $jxnQryPos =~ s/^NA://;
        my $isTopStrand = ($line[STRAND1] eq "+"); # only inline junctions handled here, both alns will be on top, or both on bottom strand
        if($lines[$jxnEdge][INSERT_SIZE] > 0){
            $lines[$jxnEdge][JXN_BASES] = 
                $isTopStrand ? 
                substr($line[SEQ1], $jxnQryPos, $lines[$jxnEdge][INSERT_SIZE]) :                 
                substr(substr($line[SEQ1], -($jxnQryPos + $lines[$jxnEdge][INSERT_SIZE])), 0, $lines[$jxnEdge][INSERT_SIZE]);
        } else {
            $lines[$jxnEdge][JXN_BASES] = "";
            my $i = 0; 
            my $alnStem = 
                $isTopStrand ? 
                substr($line[SEQ1],  $jxnQryPos) :                 
                substr($line[SEQ1], -$jxnQryPos);
            my $delBases = 
                $isTopStrand ? 
                getRefSeq($lines[$jxnEdge][CHROM1], $lines[$jxnEdge][POS1] + 1, $lines[$jxnEdge][POS2] - 1) :                 
                getRefSeq($lines[$jxnEdge][CHROM1], $lines[$jxnEdge][POS2] + 1, $lines[$jxnEdge][POS1] - 1);
            my $stemBase = substr($alnStem, $i, 1);
            while($i < $lines[$jxnEdge][EVENT_SIZE] and $stemBase eq substr($delBases, $i, 1)){
                $lines[$jxnEdge][JXN_BASES] .= $stemBase;
                $i += 1;
                $stemBase = substr($alnStem, $i, 1);
            }
            my $nJxnBases = length($lines[$jxnEdge][JXN_BASES]);
            if($nJxnBases == 0) { $lines[$jxnEdge][JXN_BASES] = "*" }            
            $lines[$jxnEdge][EVENT_SIZE] -= $nJxnBases;
            $lines[$jxnEdge][INSERT_SIZE] = -$nJxnBases;
            if($isTopStrand){
                $lines[$prevAln][POS2] += $nJxnBases; 
                $lines[$jxnEdge][POS1] = $lines[$prevAln][POS2];
            } else {
                $line[POS1] += $nJxnBases; 
                $lines[$jxnEdge][POS2] = $line[POS1];
            }
        }
        $pendingAdjustment = undef;
    }
    $prevMolId = $molId;
    push @lines, \@line;
}
$nMolecules++;
print join("", map { join("\t", @$_)."\n" } @lines);

# print summary information
printCount($nMolecules, 'nMolecules', 'input molecules');
printCount($nJunctions, 'nJunctions', 'junctions');
printCount($nInlineJunctions, 'nInlineJunctions', 'inline junctions subjected to adjustment');

1;
