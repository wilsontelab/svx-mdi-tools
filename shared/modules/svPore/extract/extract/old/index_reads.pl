use strict;
use warnings;

# reformat and index read sequences, base qualities, and moves for
#   a subset of single-alignment molecules for adapter model training
#   all SV molecules
# input is 2nd step target_reads.unaligned.bam, thus one line per a subset of reads
# sequences.txt is loaded in its entirety by analyze_edges.R and one read at a time by app

# constants
use constant {
    QNAME => 0, # sam fields
    FLAG => 1,
    RNAME => 2,
    POS => 3,
    MAPQ => 4,
    CIGAR => 5,
    RNEXT => 6,
    PNEXT => 7,
    TLEN => 8,
    SEQ => 9,
    QUAL => 10,
    TAGS => 11
};

# run the target reads to reformat (to stdout) and index (to file)
my $idxFile = "$ENV{SEQUENCES_FILE}.index";
open my $idxH, ">", $idxFile or die "could not open file: $idxFile\n";
my $offset = 0;
while(my $read = <STDIN>){
    chomp $read;
    my @read = split("\t", $read, 12);
    my ($meanBaseQual)          = ($read[TAGS] =~ m/qs:i:(\S+)/);
    my ($nSamples)              = ($read[TAGS] =~ m/ns:i:(\S+)/);
    my ($nTrimmedSamples)       = ($read[TAGS] =~ m/ts:i:(\S+)/);
    my ($channel)               = ($read[TAGS] =~ m/ch:i:(\S+)/);
    my ($pod5File)              = ($read[TAGS] =~ m/fn:Z:(\S+)/);
    my ($downsampling, $moves)  = ($read[TAGS] =~ m/mv:B:c,(\d+),(\S+)/);
    my $line = join("\t", @read[QNAME, SEQ], $moves, $downsampling, $nSamples, $nTrimmedSamples, $meanBaseQual, $channel, $pod5File)."\n"; # QUAL
    print $line;
    my $nChar = length($line);   
    print $idxH join("\t", $read[QNAME], $offset, $nChar), "\n";
    $offset += $nChar;        
}
close $idxH;

# tags in unaligned bam added by Dorado
# https://github.com/nanoporetech/bonito/blob/master/documentation/SAM.md#header
#   RG:Z:	<runid>_<basecalling_model>
#   qs:i:	mean basecall qscore rounded to the nearest integer
#   ns:i:	the number of samples in the signal prior to trimming
#   ts:i:	the number of samples trimmed from the start of the signal
#   mx:i:	read mux
#   ch:i:	read channel
#   rn:i:	read number
#   st:Z:	read start time (in UTC)
#   du:f:	duration of the read (in seconds)
#   f5:Z:	fast5 file name
#   sm:f:	scaling midpoint/mean/median (pA to ~0-mean/1-sd)
#   sd:f:	scaling dispersion (pA to ~0-mean/1-sd)
#   sv:Z:	scaling version
#   mv:B:c	sequence to signal move table (optional)
#-----------------------------------------------------------
# qs:i:18 
# du:f:45.3635 
# ns:i:181454 
# ts:i:10 
# mx:i:1 
# ch:i:2220 
# st:Z:2023-05-15T18:25:31.265+00:00 
# rn:i:26152 
# fn:Z:PAK96116_pass_d57a3f63_44212866_804.pod5 
# sm:f:104.642 
# sd:f:28.7487 
# sv:Z:quantile 
# dx:i:0 
# RG:Z:44212866912449f03e18393f276b4aed25626ba0_dna_r10.4.1_e8.2_400bps_sup@v4.1.0 
# mv:B:c,5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0 ...
