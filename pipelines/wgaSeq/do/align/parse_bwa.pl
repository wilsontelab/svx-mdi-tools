use strict;
use warnings;

# extract information that uniquely identifies each input DNA molecule
# output is one line per read pair, suitable for subsequent sorting and grouping

# initialize reporting
our $action  = "parse_BWA";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
require "$perlUtilDir/workflow.pl";
require "$perlUtilDir/sequence.pl";
require "$perlUtilDir/numeric.pl";
use vars qw($leftClip_ $rightClip_);

# environment variables
fillEnvVar(\my $nCpu,       'N_CPU');
fillEnvVar(\my $minMapQ,    'MIN_MAPQ');

# constants
use constant {
    END_READ_PAIR => 'END_READ_PAIR',
    #-------------
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
    #-------------
    ALIGN_SPLIT => 7,
    #-------------
    _PAIRED => 1,
    _PROPER_PAIR => 2,
    _UNMAPPED => 4,
    _REVERSE => 16, # SAM FLAG bits
    _SECOND => 128,
    #-------------
    SHIFT_SECOND => 7,
    #-------------
    READ1 => 0,
    READ2 => 1,
};

# working variables
my (@alns);

# process data by read pair over multiple parallel threads
launchChildThreads(\&parse_BWA);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
my $nInputReadPairs = 0;
while(my $line = <STDIN>){
    $line =~ m/^\@/ and next;
    my ($qName) = split("\t", $line, 2);     
    if($threadName and $qName ne $threadName){
        print $writeH END_READ_PAIR, "\n";
        $nInputReadPairs++;
        $writeH = $writeH[$nInputReadPairs % $nCpu + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $qName;
}
print $writeH END_READ_PAIR, "\n"; # finish last read pair
finishChildThreads();

# child process to parse BWA map output
sub parse_BWA {
    my ($childN) = @_;
    
    # auto-flush output to prevent buffering and ensure proper feed to sort
    $| = 1;

    # run BWA input one alignment at a time
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        if($line eq END_READ_PAIR){
            ($alns[READ1][0][FLAG] & _PAIRED) ? processUnmerged() : processPremerged();
            @alns = ();    
        } else { # add new alignment to growing read pair
            my @aln = (split("\t", $line, ALIGN_SPLIT))[QNAME..CIGAR];
            push @{$alns[($aln[FLAG] & _SECOND) >> SHIFT_SECOND]}, \@aln; # 0=read1, 1=read2
        }
    }
}

#---------------------------------------------------------------------
# handle read pairs that were pre-merged by fastp
#---------------------------------------------------------------------
sub processPremerged {
    
    # only consider high-quality proper molecules (same criterion as subsequent binning)
    @{$alns[READ1]} == 1 or return; # reject split reads 
    my $aln = $alns[READ1][0];        
    ($$aln[FLAG] & _UNMAPPED) and return; # reject unmapped reads (not informative)   

    # calculate position information along genome and molecule
    my $endPos = getEnd($$aln[POS], $$aln[CIGAR]);  
    my $leftClip  = $$aln[CIGAR] =~ m/$leftClip_/  ? $1 : 0;
    my $rightClip = $$aln[CIGAR] =~ m/$rightClip_/ ? $1 : 0;

    # output all fields required for grouping and binning
    print join("\t",
        $$aln[RNAME], $$aln[POS] - $leftClip, $endPos + $rightClip, # genome span covered by _entire_ SEQ
        @$aln[POS, MAPQ], $endPos - $$aln[POS] + 1 # TRUE mapped pos, MAPQ and TLEN _without_ clips
    ), "\n";   
}

#---------------------------------------------------------------------
# handle read pairs that were NOT already merged by fastp
#---------------------------------------------------------------------
sub processUnmerged {
    
    # only consider high-quality proper molecules (same criterion as subsequent binning)
    @{$alns[READ1]} == 1 or return; # reject split reads 
    @{$alns[READ2]} == 1 or return;
    my $aln1 = $alns[READ1][0];       
    my $aln2 = $alns[READ2][0];
    ($$aln1[FLAG] & _PROPER_PAIR) or return; # reject SV gaps   
    ($$aln1[FLAG] & _UNMAPPED) and return; # just in case BWA proper is iffy
    ($$aln2[FLAG] & _UNMAPPED) and return;
    $$aln1[RNAME] eq $$aln2[RNAME] or return;

    # calculate position information along genome and molecule
    ($$aln1[FLAG] & _REVERSE) and ($aln1, $aln2) = ($aln2, $aln1); # so 1st read of proper pair is Forward
    ($$aln1[FLAG] & _REVERSE) and return;
    ($$aln2[FLAG] & _REVERSE) or return;
    my $endPos = getEnd($$aln2[POS], $$aln2[CIGAR]);
    $endPos >= $$aln1[POS] or return;
    my $leftClip  = $$aln1[CIGAR] =~ m/$leftClip_/  ? $1 : 0;
    my $rightClip = $$aln2[CIGAR] =~ m/$rightClip_/ ? $1 : 0;

    # output all fields required for grouping and binning
    print join("\t",
        $$aln1[RNAME], $$aln1[POS] - $leftClip, $endPos + $rightClip, # genome span covered by _entire_ SEQ
        $$aln1[POS], min($$aln1[MAPQ], $$aln2[MAPQ]), $endPos - $$aln1[POS] + 1 # TRUE mapped pos, MAPQ and TLEN _without_ clips
    ), "\n";
}
