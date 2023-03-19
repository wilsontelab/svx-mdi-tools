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
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $EXTRACT_PREFIX,   'EXTRACT_PREFIX');
fillEnvVar(\our $ACTION_DIR,       'ACTION_DIR');
fillEnvVar(\our $N_CPU,            'N_CPU'); # user options, or derived from them
fillEnvVar(\our $WINDOW_POWER,     'WINDOW_POWER');
fillEnvVar(\our $WINDOW_SIZE,      'WINDOW_SIZE');
fillEnvVar(\our $GENOME_FASTA,     'GENOME_FASTA');
fillEnvVar(\our $USE_CHR_M,        'USE_CHR_M');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
$perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl/svLong";
map { require "$ACTION_DIR/extract/$_.pl" } qw(initialize_nodes parse_nodes);
initializeWindowCoverage();

# constants
use constant {
    END_MOLECULE => '_ERM_',
    #-------------
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
    RNAME_INDEX => 13  # added by us 
    #-------------
};

# process data by molecule over multiple parallel threads
launchChildThreads(\&parseMolecule);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
while(my $line = <STDIN>){
    $nInputAlns++;
    my ($qName) = split("\t", $line, 2);  
    if($threadName and $qName ne $threadName){
        $nInputMols++;        
        print $writeH END_MOLECULE, "\t$nInputMols\n";        
        $writeH = $writeH[$nInputMols % $N_CPU + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $qName;
}
$nInputMols++;
print $writeH END_MOLECULE, "\t$nInputMols\n";      
finishChildThreads();

# print summary information
printCount($nInputAlns, 'nInputAlns', 'input aligned segments over all molecules');
printCount($nInputMols, 'nInputMols', 'input molecules');

# child process to parse PAF molecules
sub parseMolecule {
    my ($childN) = @_;
    
    # auto-flush output to prevent buffering and ensure proper feed to sort
    $| = 1;

    # working variables
    our (@alns, 
         @nodes,    @types,    @mapQs,    @sizes,    @insSizes, 
         @alnNodes, @alnTypes, @alnMapQs, @alnSizes, @alnInsSizes, 
         $molId) = ();
  
    # run aligner output one alignment at a time
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        my @aln = split("\t", $line, 13); 

        # process all alignments from one source molecule    
        if($aln[0] eq END_MOLECULE){

            # characterize the path of a single molecule
            $molId = ($aln[1] * 100) + $childN;
            @alns = sort { $$a[QSTART] <=> $$b[QSTART] } @alns; # sort alignments in query order (could be right to left on bottom strand)

            if(@alns == 2){

            foreach my $i(0..$#alns){

                $i > 0 and printSplitJunction($alns[$i-1], $alns[$i]);

                # my ($cigar) = ($alns[0][PAF_TAGS] =~ m/cg:Z:(\S+)/);
                # $cigar =~ m/\d+D\d+I/ or next;

                (@alnNodes, @alnTypes, @alnMapQs, @alnSizes, @alnInsSizes) = ();
                processAlignedSegment($alns[$i]);
                push @nodes,    @alnNodes;
                push @types,    @alnTypes;
                push @mapQs,    @alnMapQs;
                push @sizes,    @alnSizes;
                push @insSizes, @alnInsSizes;
            }
            
            
            if(@nodes){

            # set junction MAPQ as minimum MAPQ of flanking alignments
            if(@mapQs > 1){
                for (my $i = 1; $i <= $#mapQs - 1; $i += 2){
                    $mapQs[$i] = min($mapQs[$i - 1], $mapQs[$i + 1])
                }
            }

            # print one line per node pair in the collapsed molecule sequence
            # print en bloc to keep together when merging parallel streams
            my $out = "";  
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
            print $out, "\n";    

            }
            }

            # reset for next molecule
            (@alns, @nodes, @types, @mapQs, @sizes, @insSizes) = ();

        # add new alignment to growing source molecule    
        } else{
            $aln[RNAME_INDEX] = $chromIndex{$aln[RNAME]} || 0; # unknown sequences places into chrom/window zero 
            push @alns, \@aln;
        }
    }

    # print the partial coverage map from this thread
    printWindowCoverage($childN);
}
