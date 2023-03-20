use strict;
use warnings;

# examine all SV junctions for evidence of chimeric reads
# see:
#   https://github.com/nanoporetech/duplex-tools/blob/master/fillet.md
#   https://github.com/nanoporetech/duplex-tools/blob/master/duplex_tools/split_on_adapter.py
# occurs when a molecule has duplex reads, or when two molecules pass one after another
# when this occurs, the read structure is 
#   >>====<<>>=====<<
# where >> and << are the adapter and its complement
# Guppy (I think) trims outer adapters, but not internal adapters.
# duplex_tools can split on internal adapters, but prefer to do it here, with alignment guidance
# svLong wants to collapse, not keep, redundant strands in duplex pairs
# note that this procedure applies to CIGAR indels also, not just split alignments

# constants
use constant {
    QNAME => 0, # PAF fields
    QLEN => 1,
    QSTART => 2,
    QEND => 3,
    STRAND => 4,
    RNAME => 5,
    RLEN => 6,
    RSTART => 7,
    REND => 8,
    N_MATCHES => 9,
    N_BASES => 10,
    MAPQ => 11,
    PAF_TAGS => 12,
    RNAME_INDEX => 13,  # added by us  
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    #-------------
    MATCH_OP      => "M",
    NULL_OP       => "X"
};

# working variables
use vars qw($INPUT_DIR $molId @nodes @types @mapQs @sizes @insSizes @outAlns);   

# my $mask_size_default_head = 5;
# my $mask_size_default_tail = 14;
# my $mask_size_default_N = 11;
my $HEAD_ADAPTER = 'AATGTACTTCGTTCAGTTACGTATTGCT';
my $TAIL_ADAPTER = 'GCAATACGTAACTGAACGAAGT';
my $adapterTarget = $TAIL_ADAPTER.$HEAD_ADAPTER;
# my $rcAdapterTarget = rc($adapterTarget);

sub checkJunctionsForAdapters {
    my $nJxns = @types;

    # check whether the molecule is consistent with a duplex foldback inversion
    # if yes, keep one half only and mark it as duplex
    # self-symmetry is key to making this assessment
    my @inversionIs;
    map { $types[$_] eq INVERSION and push @inversionIs, $_ } 0..$#types;
    if(@inversionIs % 2 == 1){ # demand an odd number of inversions, so the fold-back inversion can be central
        printMolecule();
        my @read = getIndexedRead($molId);

        my $refPos11 = $outAlns[0][STRAND] eq "+" ? $outAlns[0][RSTART] + 1 : $outAlns[0][REND];
        my $refPos12 = $outAlns[0][STRAND] eq "+" ? $outAlns[0][REND]       : $outAlns[0][RSTART] + 1;
        my $refPos21 = $outAlns[2][STRAND] eq "-" ? $outAlns[2][REND]       : $outAlns[2][RSTART] + 1;        
        my $refPos22 = $outAlns[2][STRAND] eq "-" ? $outAlns[2][RSTART] + 1 : $outAlns[2][REND];

        my $qryPos11 = $outAlns[0][QSTART];
        my $qryPos12 = $outAlns[0][QEND];
        my $qryPos21 = $outAlns[2][QSTART];        
        my $qryPos22 = $outAlns[2][QEND]; 

        print join("\t", $qryPos11, $qryPos12, $qryPos21, $qryPos22), "\n";
        print join("\t", $refPos11, $refPos12, $refPos21, $refPos22), "\n";
        $read[0] =~ m/^(\S+)/;
        print "---".$1."===\n";
        print "---".length($read[1])."===\n";        
        print "---".$read[2]."===\n";
        print "---".length($read[3])."===\n"; 
        print join("\t", @{$outAlns[0]}[QNAME..REND]), "\n";

        print substr($read[1], $qryPos12 - 20, 100), "\n";

        # 
        # foreach my $i(@inversionIs){
            
        # }
    }

    # check each junction for the potential that it is a chimera
    # if yes, split the molecule into two
    # does it have an expected insertion size?
    # does that insertion match $adapterTarget? (use indexed retrieval)
}

# print a molecule to STDOUT
sub printMolecule {
    my $out = "----------------------------\n";  
    foreach my $i(1..$#nodes){
        $out .= join(
            "\t", 
            $molId,
            $nodes[$i-1], 
            $nodes[$i], 
            $types[$i-1], 
            $mapQs[$i-1], 
            $sizes[$i-1], 
            $insSizes[$i-1]
        )."\n";
    }
    print $out; #, "\n"; 
}

1;
