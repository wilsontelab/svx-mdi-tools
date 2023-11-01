use strict;
use warnings;

# load dependencies
our $script = "extract_adapters";
our $error  = "$script error";
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);

# constants
use constant {
    END_READ_EDGES => "__END_READ_EDGES__",
    #=============
    QNAME => 0, # edge format fields
    NODE1 => 1,
    QSTART => 2,    
    NODE2 => 3,
    QEND => 4,
    MAPQ => 5,
    CIGAR => 6,
    GAP_COMPRESSED_IDENTITY => 7,
    EDGE_TYPE => 8,
    EVENT_SIZE => 9,
    INSERT_SIZE => 10,
    N_STRANDS => 11,
    #=============
    _QNAME => 0, # SAM fields
    _FLAG => 1,
    _RNAME => 2,
    _POS => 3,
    _MAPQ => 4,
    _CIGAR => 5,
    _RNEXT => 6,
    _PNEXT => 7,
    _TLEN => 8,
    _SEQ => 9,
    _QUAL => 10,
    _TAGS => 11,
    #=============
    ALIGNMENT     => "A"
};

# environment variables
fillEnvVar(\our $ONT_LIBRARY_TYPE,  'ONT_LIBRARY_TYPE');
fillEnvVar(\our $UBAM_DIR,          'UBAM_DIR');

# after scoping out the range of true adapter lengths, set them here for sequence exploration
my $minClip5 = $ONT_LIBRARY_TYPE eq "ligation" ? 27 : 104;
my $maxClip5 = $ONT_LIBRARY_TYPE eq "ligation" ? 41 : 128;
my $minClip3 = $ONT_LIBRARY_TYPE eq "ligation" ? 9  : 1e8; # there is no 3' adapter for rapid kit
my $maxClip3 = $ONT_LIBRARY_TYPE eq "ligation" ? 15 : 0;
my $maxLength5 = 41; # only report the part of the adapter attached to gDNA at its 3' end

# open the unaligned read files for parallel threading of QNAMES
open my $ubamH, "-|", "samtools cat $UBAM_DIR/*.unaligned.bam | samtools view -" or die "could not open bam files in $UBAM_DIR: $!\n";

# pull putative candidate adapters on the ends of simple alignments
my ($prevQName, @edgeArrays) = (0);
sub printAdapters {
    my $isSv = (@edgeArrays > 1);
    my $read = <$ubamH>;
    while($read !~ m/^$prevQName/){ # handle reads found in SAM that failed to align and are absent from PAF
        $read = <$ubamH>;
        $read or last;
    }
    if(!$read){
        die "$error:\ngot to end of reads in:\n$UBAM_DIR\nwithout encountering read ID:\n$prevQName\n";
    }
    $isSv and return;
    my @read = split("\t", $read, 12);
    $read[_TAGS] =~ m/dx:i:(\S+)/;
    my $dx = $1;
    $dx == 0 or return; # see below, duplex reads often don't have adapters, simplex more reliable
    my $readLen = length($read[_SEQ]);
    my $adapterLen5 = $edgeArrays[0][QSTART];
    my $adapterLen3 = $readLen - $edgeArrays[0][QEND];
    my $hasAdapter5 = ($adapterLen5 >= $minClip5 and $adapterLen5 <= $maxClip5);
    my $hasAdapter3 = ($adapterLen3 >= $minClip3 and $adapterLen3 <= $maxClip3);
    print join("\t", 
        $adapterLen5,
        $adapterLen3,
        $hasAdapter5 ? substr((" " x ($maxClip5 - $adapterLen5)).substr($read[_SEQ], 0, $adapterLen5), -$maxLength5) : "NA",
        $hasAdapter3 ? substr($read[_SEQ], -$adapterLen3)   : "NA"
    ), "\n";
}
while(my $edge = <STDIN>){
    chomp $edge;
    my @edgeArray = split("\t", $edge);
    if($prevQName and $prevQName ne $edgeArray[QNAME]){
        printAdapters();
        @edgeArrays = ();
    }
    push @edgeArrays, \@edgeArray;
    $prevQName = $edgeArray[QNAME];
}
close $ubamH;

# ==================================
# adapter lengths (distribution tails trimmed)
# ==================================
# ONT Ligation kit
# ==================================
# 5' end adapter lengths, SIMPLEX ONLY
# ----------------------------------
# 0       1853 **** small minority of reads with no adapter, likely capture randomly by a motor?
# 1       151
# 2       106
# 3       110
# ...
# 21      125
# 22      207
# 23      276
# 24      370
# 25      608
# 26      995
# 27      1400
# 28      1674
# 29      2236
# 30      3009
# 31      4396
# 32      5689 ****, the requisite duplex stem of the adapter plus the motor loading sequence
# 33      5358       << in net, easy to train and find these adapters >> 
# 34      5150
# 35      4596
# 36      4671
# 37      4322
# 38      3423
# 39      2395
# 40      1631
# 41      1118
# ----------------------------------
# 3' end adapter lengths, SIMPLEX ONLY
# ----------------------------------
# 0       5144 ****, a smaller subset that runs out of the gDNA with no adapter
# 1       3702
# 2       4460
# 3       2385
# 4       1430
# 5       1024
# 6       847
# 7       806
# 8       816
# 9       1097
# 10      2157
# 11      4485
# 12      9204 ****, a larger subset with a short stem from the duplex portion of the adapter
# 13      4685       << in net, can train and find these adapters stems, but expect higher failure rate >>
# 14      3597
# 15      1411
# 16      882
# 17      606
# 18      366
# 19      257
# 20      232
# 21      212
# 22      202
# 23      172
# 24      140
# 25      109
# 26      126
# 27      129
# ----------------------------------
# 5' end adapter lengths, DUPLEX ONLY
# ----------------------------------
# 0       2241 ****, a larger subset where the duplex sequence is incomplete and therefore lacks an adapter
# 1       108
# 2       100
# ...
# 24      102
# 25      107
# 26      149
# 27      158
# 28      173
# 29      198 ****, a smaller subset with a complete duplex sequence that includes the 5' adapter
# 30      160       << in net, most (but not all) ligation kit duplex ends lack adapters >>
# 31      169
# 32      172
# 33      140
# 34      109
# 35      82
# 36      71
# 37      52
# 38      34
# 39      20
# 40      11
# ----------------------------------
# 3' end adapter lengths, DUPLEX ONLY
# ----------------------------------
# 0       4751 ****, essentially all reads run out of the duplex at the 3' end with no adapter
# 1       241
# 2       32
# 3       36
# 4       16
# 5       10
# ----------------------------------
# most frequent 5' adapters
# ----------------------------------
    #      ATGTCCTGT ACTTCGTTCAGTTACGTATTGCT       300
    #       TGTCCTGT ACTTCGTTCAGTTACGTATTGCT       108
    #     TATGTCCTGT ACTTCGTTCAGTTACGTATTGCT       96
    #           ATGT ACTTCGTTCAGTTACGTATTGCT       93
    #    TTATGTCCTGT ACTTCGTTCAGTTACGTATTGCT       90
# ==================================
# ONT Rapid kit
# ==================================
# 5' end adapter lengths, SIMPLEX ONLY
# ----------------------------------
# 0       67
# 1       2
# ...
# 104     1101
# 105     1426
# 106     1813
# 107     2178
# 108     2670
# 109     3174
# 110     3751
# 111     4182
# 112     4638
# 113     5041
# 114     5414
# 115     5642
# 116     5708 *****, the requisite composite adapter, expected to include the mosaic end bonded to gDNA plus barcodes, etc.
# 117     5613        << main point of training is to located the junction to the mosaics ends >>
# 118     5328
# 119     4944
# 120     4564
# 121     4046
# 122     3555
# 123     2904
# 124     2357
# 125     1943
# 126     1640
# 127     1246
# 128     1007
# ----------------------------------
# 3' end adapter lengths, SIMPLEX ONLY
# ----------------------------------
# 0       27744 *****, 3' ends have no adapter, since mosaic end not bonded to gDNA by Tn5
# 1       19150        << therefore, skip adapter training and finding on rapid kit 3' ends >>
# 2       18205
# 3       9663
# 4       5614
# 5       3606
# 6       2538
# 7       1817
# 8       1241
# ----------------------------------
# 5' end adapter lengths, DUPLEX ONLY
# ----------------------------------
# 0       91 ****, ~all reads run out of the duplex at the 5' end with no adapter, since the corresponding 3' end had no adapter
# 1       3        << in net, duplex reads with the rapid kit ~never have adapter sequences >>
# 3       3
# 4       3
# 5       1
# 6       1
# 8       1
# 9       1
# 10      1
# 17      1
# 103     1
# 106     1
# 118     1
# ----------------------------------
# 3' end adapter lengths, DUPLEX ONLY
# ----------------------------------
# 0       100 ****, ~all reads run out of the duplex at the 3' end with no adapter, c/w above
# 1       5
# 2       1
# 3       2
# 5       1
# ----------------------------------
# most frequent 5' adapters
# ----------------------------------
#     TTTATCGTGAAACGCT TTCGCGTTTTTCGTGCGCCGCTTCA       21099
#    ATTTATCGTGAAACGCT TTCGCGTTTTTCGTGCGCCGCTTC       6304
#   CATTTATCGTGAAACGCT TTCGCGTTTTTCGTGCGCCGCTT       1724
#    ATTTATCGTGAAACGCT TTCGCGTTTTTCGTGCGCCGCT-CA       1634
# CGCATTTATCGTGAAACGCT TTCGCGTTTTTCGTGCGCCGC       1198
# ==================================
