use strict;
use warnings;

# execute the multi-step process of consensus making
#   STEP 1 - establish coherent merge levels of all molecules in the group
#   STEP 2 - if needed, downsample read pairs per molecule+strand to a reasonable number
#   STEP 3 - make a single-strand consensus for each read, or merged read, of each source strand
#   STEP 4 - make a duplex consensus from the single-strand consensuses when available
# a list of QNAMEs is saved for recovery from primary bam
# all merged or consensus molecules are passed forward for genome remapping
# genomic bases from source molecules should mostly be represented once in output reads
# non-genomic bases (UMIs, adapters) should be clipped away in advance of next mapping
# remaining large clips should mainly correspond to SVs at:
#   the outer molecule end
#   the inner side of read pairs that do not overlap
    
# initialize reporting
our $script  = "make_consensus";
our ($nReadPairs, $nMolecules,
     $nThreadMolecules, $nThreadMerged, $nThreadDuplex) = (0) x 100;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general IUPAC smith_waterman);
resetCountFile();

# environment variables
fillEnvVar(\my $N_CPU,                'N_CPU');
fillEnvVar(\my $ACTION_DIR,           'ACTION_DIR');
fillEnvVar(\my $CONSENSUS_PREFIX,     'CONSENSUS_PREFIX');
fillEnvVar(\my $DOWNSAMPLE_N,         'DOWNSAMPLE_N');
fillEnvVar(\our $MIN_SW_SCORE_FACTOR, 'MIN_SW_SCORE_FACTOR');
fillEnvVar(\our $CONSENSUS_FACTOR,    'CONSENSUS_FACTOR');

# load additional dependencies
map { require "$ACTION_DIR/consensus/$_.pl" } qw(merge consensus quality);

# constants
use constant {
    READ_PAIR_ID => 0, # columns in read-pair lines
    MOL_STRAND => 1,
    SEQ1 => 2, # SEQ indices refer to reads
    QUAL1 => 3,
    SEQ2 => 4,
    QUAL2 => 5,
    QNAME => 6,
    SEQ_MERGED => 7, # buffers used by consensus scripts
    QUAL_MERGED => 8,
    #-----------------
    MOL_MARKER => 0, # columns in molecule lines
    MOL_ID => 1,
    UMI1 => 2,
    UMI2 => 3,    
    STRAND_COUNT1 => 4, # indices refer to strands
    STRAND_COUNT2 => 5,
    IS_DUPLEX => 6, # buffers used by consensus scripts
    IS_MERGED => 7,
    #-------------------
    STRAND1 => 0, # consensus array indices, for code readability
    STRAND2 => 1, # so can have combinations like STRAND1|READ1|SEQ or DUPLEX|MERGED|QUAL
    DUPLEX => 2,    
    READ1 => 0,
    READ2 => 1,
    MERGED => 2,
    SEQ => 0,
    QUAL => 1,
};

# loop the input, breaking into chunks corresponding to individual source molecules
# process data by molecule over multiple parallel threads
launchChildThreads(\&makeConsensuses);
use vars qw(@readH @writeH);
my $writeH = $writeH[$nMolecules % $N_CPU + 1];
while(my $line = <STDIN>){
    print $writeH $line;
    if($line =~ m/^\@M/){
        $nMolecules++;
        $writeH = $writeH[$nMolecules % $N_CPU + 1];
    } else {
        $nReadPairs++;
    }
}
finishChildThreads();

# print summary information
printCount($nReadPairs,     'nReadPairs', 'input read pairs');
printCount($nReadPairs * 2, 'nReads',     'input reads (nReadPairs * 2)');
printCount($nMolecules,     'nMolecules', 'input source DNA molecules');

# child process to parse read-pairs into a (merged) consensus
sub makeConsensuses {
    my ($childN) = @_;

    # working variables
    our (@readPairs, $mol, @strands, @downsample, @consensus) = ();
    
    # output file handles
    my $fqFile = "$CONSENSUS_PREFIX.$childN.fq.gz";
    open my $fqH, "|-", "gzip -c > $fqFile" or die "could not open $fqFile: $!\n";
    my $nameMapFile = "$CONSENSUS_PREFIX.$childN.name_map.gz";
    open my $nameMapH, "|-", "gzip -c > $nameMapFile" or die "could not open $nameMapFile: $!\n";

    # run the read pairs
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        my @line = split("\t", $line);
    
        # molecule line, finishes a set of read pairs
        if($line[0] eq '@M'){ 
            $mol = \@line;
            $nThreadMolecules++;

            # STEP 1 - establish coherent merge levels of all molecules in the group
            @strands = map { $$mol[$_] ? $_ - STRAND_COUNT1 : () } STRAND_COUNT1, STRAND_COUNT2;              
            parseMergeLevels(); # may possibly override @strands
            $$mol[IS_MERGED] and $nThreadMerged++;
            my @outReads = $$mol[IS_MERGED] ? (MERGED) : (READ1, READ2);

            # determine needed strandedness information
            $$mol[IS_DUPLEX] = @strands - 1;
            $$mol[IS_DUPLEX] and $nThreadDuplex++;
            my $molName = join(":",
                @$mol[MOL_ID, UMI1, UMI2, IS_MERGED, # in parallel order as genomex-mdi-tools align, but IS_MERGED = (0,1)
                      IS_DUPLEX, STRAND_COUNT1, STRAND_COUNT2], # new values appended to QNAME
            );            

            # STEP 2 - if needed, downsample read pairs per molecule+strand to a managable but informative number
            foreach my $strand(@strands){              
                @{$downsample[$strand]} = @{$readPairs[$strand]} > $DOWNSAMPLE_N ?
                    0..($DOWNSAMPLE_N - 1) : # group_reads.pl randomized read pairs within the molecule group
                    0..$#{$readPairs[$strand]};
            }
  
            # STEP 3 - make a single-strand consensus for each read, or merged read, of each source strand
            @consensus = ();
            foreach my $strand(@strands){
                if($$mol[IS_MERGED]){
                    makeStrandMergeConsensus($strand);
                } else {
                    makeStrandReadConsensus($strand);
                }   
            }
            
            # STEP 4 - make a duplex consensus from the single-strand consensuses when available
            $$mol[IS_DUPLEX] and makeDuplexConsensus(@outReads);
            
            # STEP LAST - print the final consensus information in FASTQ format for remapping
            # print either 1 or 2 interleaved reads per molecule depending on whether read-pair merging occurred
            # those reads both come from either one strand (simplex) or one duplex consensus
            # and represent either one merged read, two unmerged reads, or an orpan unmerged whose partner failed consensus
            # duplex reads that failed duplex consensus building are not printed
            my $outStrand = $$mol[IS_DUPLEX] ? DUPLEX : $strands[0];    
            foreach my $read(@outReads){ 
                if($consensus[$outStrand] and
                   $consensus[$outStrand][$read] and
                   $consensus[$outStrand][$read][SEQ]){
                    print $fqH join("\n",
                        '@'.$molName,
                        $consensus[$outStrand][$read][SEQ],
                        "+",
                        $consensus[$outStrand][$read][QUAL]
                    ), "\n";
                }  
            }
            
            # save a table that correlates all primary bam QNAMEs to consensus molecule names
            # format: molName QNAME1,QNAME2,QNAME3[,...]
            print $nameMapH join("\t", $molName, join(",", map {
                map { $$_[QNAME] } @{$readPairs[$_]}
            } @strands)), "\n";
            
            # prepare for the next molecule
            @readPairs = ();
            
        # read pair line, with one claimed read pair of a molecule+strand
        # add to the growing molecule+strand group
        } else {
            push @{$readPairs[$line[MOL_STRAND]]}, \@line;
        }
    }
    
    # finish thread
    close $fqH;
    close $nameMapH;
    
    # print molecule counts
    printCount($nThreadMolecules, "nThreadMolecules-$childN", "total molecules processed in thread $childN");
    printCount($nThreadMerged,    "nThreadMerged-$childN",    "molecules were merged in thread $childN");
    printCount($nThreadDuplex,    "nThreadDuplex-$childN",    "molecules were duplex in thread $childN");
}
