use strict;
use warnings;

# extract SV information from collated, name-sorted long-read alignments

# initialize reporting
our $script = "extract_nodes";
our $error  = "$script error";
my ($nInputAlns, $nInputMols) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
our ($matchScore, $mismatchPenalty, $gapOpenPenalty, $gapExtensionPenalty) = 
    (1,           -1.5,             -1.5,            -2);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman); # faidx 
resetCountFile();

# environment variables
fillEnvVar(\our $EXTRACT_PREFIX,   'EXTRACT_PREFIX');
fillEnvVar(\our $PIPELINE_DIR,     'PIPELINE_DIR');
fillEnvVar(\our $GENOMEX_MODULES_DIR, 'GENOMEX_MODULES_DIR');
fillEnvVar(\our $WINDOW_SIZE,      'WINDOW_SIZE');
fillEnvVar(\our $MIN_SV_SIZE,      'MIN_SV_SIZE');
fillEnvVar(\our $GENOME_FASTA,     'GENOME_FASTA');
our $USE_CHR_M = 1;

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
# require "$GENOMEX_MODULES_DIR/align/dna-long-read/get_indexed_reads.pl";
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svPore";
map { require "$PIPELINE_DIR/extract/$_.pl" } qw(initialize_windows parse_nodes);
$perlUtilDir = "$ENV{MODULES_DIR}/parse_nodes";
map { require "$perlUtilDir/$_.pl" } qw(parse_nodes_support);
initializeWindowCoverage();

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
};

# working variables
our (@alns, $molId,
     @nodes,    @mapQs,    @cigars,    @alnQs,    @types,    @sizes,    @insSizes,    @outAlns,
     @alnNodes, @alnMapQs, @alnCigars, @alnAlnQs, @alnTypes, @alnSizes, @alnInsSizes, @alnAlns) = ();

# process data by molecule over multiple parallel threads
my ($prevQName);
$| = 1;
while(my $line = <STDIN>){
    $nInputAlns++;
    my @aln = split("\t", $line, 13);
    $aln[RNAME_INDEX] = $chromIndex{$aln[RNAME]} || 0; # unknown sequences places into chrom/window zero 
    if($prevQName and $aln[QNAME] ne $prevQName){
        parseMolecule();
    }    
    push @alns, \@aln;
    $prevQName = $aln[QNAME];
}
parseMolecule(); 

# print summary information
printCount($nInputMols, 'nInputMols', 'input molecules');
printCount($nInputAlns, 'nInputAlns', 'input aligned segments over all molecules');

# child process to parse PAF molecules
sub parseMolecule {
    $nInputMols++;
    $molId = $alns[0][QNAME];   

    # characterize the path of a single molecule
    @alns = sort { $$a[QSTART] <=> $$b[QSTART] } @alns; # sort alignments in query order (could be right to left on bottom strand)
    foreach my $i(0..$#alns){

        # add information on the junction between two alignments
        if($i > 0){
            my $jxn = processSplitJunction($alns[$i-1], $alns[$i]);
            push @mapQs,    0;
            push @cigars,   "NA";
            push @alnQs,    0;
            push @types,    $$jxn{jxnType};            
            push @sizes,    $$jxn{svSize};
            push @insSizes, $$jxn{insSize};
            push @outAlns,  [];
        }         
        # add each alignment    
        (@alnNodes, @alnMapQs, @alnCigars, @alnAlnQs, @alnTypes, @alnSizes, @alnInsSizes, @alnAlns) = ();
        processAlignedSegment($alns[$i]);
        push @nodes,    @alnNodes;
        push @mapQs,    @alnMapQs;
        push @cigars,   @alnCigars;
        push @alnQs,    @alnAlnQs;        
        push @types,    @alnTypes;        
        push @sizes,    @alnSizes;
        push @insSizes, @alnInsSizes;
        # map { 
        #     if($alnInsSizes[$_] !~ m/\t/){ # add query positions in xStart and xEnd in CIGAR junctions (hopefully will be few of these)
        #         my $nodePos1 = $alnAlns[$_ - 1]         eq "+" ? $alnAlns[$_ - 1][REND] : $alnAlns[$_ - 1][RSTART] + 1;
        #         my $nodePos2 = $alnAlns[$_ + 1][STRAND] eq "-" ? $alnAlns[$_ + 1][REND] : $alnAlns[$_ + 1][RSTART] + 1;
        #         $alnInsSizes[$_] = join("\t", $alnInsSizes[$_], $nodePos1, $nodePos2, FROM_CIGAR);
        #     }
        #     $alnInsSizes[$_];
        # } 0..$#alnInsSizes;
        push @outAlns,  @alnAlns;
    }

    # examine SV junctions for evidence of duplex reads
    # adjusts the output arrays as needed
    my $nStrands = checkForDuplex(scalar(@types));

    # set junction MAPQ as minimum MAPQ of the two flanking alignments
    fillJxnQs();

    # print one line per node pair, i.e., per edge, in the collapsed molecule sequence
    printMolecule($molId, $nStrands);

    # reset for next molecule
    (@alns, 
     @nodes,    @mapQs,    @cigars,    @alnQs,    @types,    @sizes,    @insSizes,    @outAlns) = ();
}
