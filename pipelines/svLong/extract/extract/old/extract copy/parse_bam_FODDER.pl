use strict;
use warnings;

# extract information that summarizes the nature of amplicon molecules

# initialize reporting
our $script = "parse_bam";
my ($nInputAlns, $nReadPairs) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\my $N_CPU,        'N_CPU');

# constants
use constant {
    END_READ_PAIR => '_ERP_',
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
    # #-------------
    REVERSE => 11,  # additional fields to take working values
    STRAND => 12,    
    CLIP1 => 13, # here, 1/2 suffix refers to the read number the alignment is "acting as"
    CLIP2 => 14,
    GROUP_POS1 => 15,
    GROUP_POS2 => 16,
    SIDE1 => 17,
    SIDE2 => 18,
    #-------------
    ALIGN_SPLIT => 12,
    #-------------
    _REVERSE => 16,
    #-------------
    SHIFT_REVERSE => 4, # how far to shift those bits to yield binary values
    SHIFT_SECOND => 7,
    #-------------
    READ1 => 0,
    READ2 => 1,
};

# regular expressions for UMI and clip handling
use vars qw($leftClip_ $rightClip_);
my %actingIs = (CLIP0 => CLIP1, # so, e.g., "CLIP".READ1 => "CLIP1"
                CLIP1 => CLIP2,
                GROUP_POS0 => GROUP_POS1,
                GROUP_POS1 => GROUP_POS2,
                SIDE0 => SIDE1,
                SIDE1 => SIDE2);

# working variables
# my $nullSymbol = '*';
my @topStrandIs = (READ1, READ2, GROUP_POS1, GROUP_POS2, CLIP1, CLIP2, SIDE1, SIDE2);
my @botStrandIs = (READ2, READ1, GROUP_POS2, GROUP_POS1, CLIP2, CLIP1, SIDE2, SIDE1);
my (@alns, @umis, @refAlns, @end1Alns) = ();

# process data by read pair over multiple parallel threads
launchChildThreads(\&parseReadPair);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
while(my $line = <STDIN>){
    $nInputAlns++;
    my ($qName) = split("\t", $line, 2);  
    # name = molId:ampliconId:nOverlapBases:molCount:merged:readN 
    $qName =~ m/(.+):\d/; # strip the trailing readN to group read by readPair
    my $pairName = $1;
    if($threadName and $pairName ne $threadName){   
        $nReadPairs++;     
        print $writeH END_READ_PAIR, "\n";
        $writeH = $writeH[$nReadPairs % $N_CPU + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $pairName;
}
$nReadPairs++; 
print $writeH END_READ_PAIR, "\n"; # finish last read pair
finishChildThreads();

# print summary information
printCount($nInputAlns, 'nInputAlns', 'input aligned segments over all read pairs');
printCount($nReadPairs, 'nReadPairs', 'input read pairs');

# child process to parse bam read pairs
sub parseReadPair {
    my ($childN) = @_;
    
    # auto-flush output to prevent buffering and ensure proper feed to sort
    $| = 1;

    # run aligner output one alignment at a time
    my $readH = $readH[$childN];
    my ($molId, $ampliconId, $nOverlapBases, $molCount, $isMerged);
    while(my $line = <$readH>){
        chomp $line;
        if($line eq END_READ_PAIR){
            
            # extract information on the source read pair
            # name = molId:ampliconId:nOverlapBases:molCount:merged:readN
            ($molId, $ampliconId, $nOverlapBases, $molCount, $isMerged) = $alns[READ1] ? split(":", $alns[READ1][0][QNAME]) : ();

            # discard orphan reads since cannot associate UMIs with endpoints
            # process others according to current merge state
            if($isMerged or ($alns[READ1] and $alns[READ2])){ 
                $isMerged ? 
                    processPremerged($molId, $ampliconId, $nOverlapBases, $molCount) : 
                    processUnmerged( $molId, $ampliconId, $nOverlapBases, $molCount);
            }

            # prepare for next read-pair
            @alns = ();    
            
        } else { # add new alignment to growing read pair
            my @aln = (split("\t", $line, ALIGN_SPLIT))[QNAME..QUAL];
            $aln[QNAME] =~ m/(.+):(\d)/; # strip the trailing readN to group read by readPair
            $aln[QNAME] = $1;
            push @{$alns[ $2 - 1 ]}, \@aln; # 0=read1, 1=read2
        }
    }
}

#---------------------------------------------------------------------
# handle read pairs that were pre-merged by fastp
#---------------------------------------------------------------------
sub processPremerged {
    my ($molId, $ampliconId, $nOverlapBases, $molCount) = @_;
    @refAlns = ();

    # calculate clip lengths for all alignments
    map {
        my $clip1  = ($$_[FLAG] & _REVERSE) ? $rightClip_ : $leftClip_;
        my $clip2  = ($$_[FLAG] & _REVERSE) ? $leftClip_ : $rightClip_;
        $$_[CLIP1] = $$_[CIGAR] =~ m/$clip1/ ? $1 : 0;
        $$_[CLIP2] = $$_[CIGAR] =~ m/$clip2/ ? $1 : 0;
    } @{$alns[READ1]};    

    # most reads have only one alignment, use it as reference alignment at both ends
    if(@{$alns[READ1]} == 1){
        @refAlns = ($alns[READ1][0], $alns[READ1][0]);    
 
    # otherwise, sort the alignments in order from read 1 to read 2 along the read pair
    # shortest clip1 is closest to READ1 primer
    } else {
        my @sorted = sort {$$a[CLIP1] <=> $$b[CLIP1]} @{$alns[READ1]};
        $refAlns[READ1] = $sorted[0]; # amplicons always use the outermost alignment as reference, regardless of MAPQ    
        $refAlns[READ2] = $sorted[$#sorted];  
    }

    # set strandedness information
    setAlignmentStrands();

    # set hypothetical outermost mapped position as used for read grouping
    setMergedGroupPos($refAlns[READ1], READ1);
    setMergedGroupPos($refAlns[READ2], READ2);
    
    # order duplex endpoints to set molecule strand
    my ($read1, $read2, $groupPos1, $groupPos2, $clip1, $clip2, $side1, $side2) = sortReferenceAlignments();

    # # determine which alignment yields the reference genome sequence OUTSIDE of inverted segments
    # # remember that aligner has already RC'ed REVERSE alignments
    # my $seqRead = (
    #     $refAlns[$read1][STRAND] == $refAlns[$read2][STRAND] or
    #     ($read1 == READ1 and !$refAlns[READ1][STRAND]) or
    #     ($read1 == READ2 and  $refAlns[READ2][STRAND])
    # ) ? $read1 : $read2;

    # output the molecule summary
    print summarizeMolecule(
        $molId, $ampliconId, $nOverlapBases, $molCount, 1,
        $read1, $read2, $side1, $side2, $groupPos1, $groupPos2
    ), "\n";  
}
sub setMergedGroupPos {
    my ($aln, $actingRead) = @_;
    my $groupPosI = $actingIs{"GROUP_POS".$actingRead};
    my $clipI     = $actingIs{"CLIP".$actingRead};
    my $sideI     = $actingIs{"SIDE".$actingRead};
    if(($actingRead == READ1 and !$$aln[REVERSE]) or
       ($actingRead == READ2 and  $$aln[REVERSE])){
        $$aln[$groupPosI] = $$aln[POS] - $$aln[$clipI];
        $$aln[$sideI] = "R";
    } else {
        $$aln[$groupPosI] = getEnd($$aln[POS], $$aln[CIGAR]) + $$aln[$clipI];
        $$aln[$sideI] = "L";
    } 
}

#---------------------------------------------------------------------
# handle read pairs that were NOT already merged by fastp
# at entry, could be either pre-merge failures or unmergable (non-overlapping) pairs
#---------------------------------------------------------------------
sub processUnmerged {
    my ($molId, $ampliconId, $nOverlapBases, $molCount) = @_;
    (@refAlns, @end1Alns) = ();    

    # analyze each read of the pair
    foreach my $read(READ1, READ2){

        # calculate outer clip lengths on all alignments
        # use to sort alignment and to determine grouping position
        my $clipI = $actingIs{"CLIP".$read};          
        map {
            my $clip = ($$_[FLAG] & _REVERSE) ? $rightClip_ : $leftClip_;
            $$_[$clipI] = $$_[CIGAR] =~ m/$clip/ ? $1 : 0;
        } @{$alns[$read]};

        # most reads have only one alignment, use it
        if(@{$alns[$read]} == 1){
            push @refAlns,  $alns[$read][0];
            push @end1Alns, $alns[$read][0];
        
        # otherwise, sort the alignment in molecule order
        # shortest outer clip is outermost alignment on the read
        } else {           
            my @sorted = sort {$$a[$clipI] <=> $$b[$clipI]} @{$alns[$read]};
            push @end1Alns, $sorted[$read==READ1 ? 0 : $#sorted];                 
            push @refAlns, $sorted[0]; # amplicons always use the outermost alignment as reference, regardless of MAPQ    
        } 
    }
    
    # set strandedness information; REVERSE as aligned, STRAND corrected to source molecule
    setAlignmentStrands();
    $refAlns[READ2][STRAND] = ($refAlns[READ2][STRAND] + 1) % 2;

    # set hypothetical outermost mapped position as used for read grouping
    setUnmergedGroupPos($refAlns[READ1], READ1);
    setUnmergedGroupPos($refAlns[READ2], READ2);

    # order duplex endpoints to set molecule strand
    my ($read1, $read2, $groupPos1, $groupPos2, $clip1, $clip2, $side1, $side2) = sortReferenceAlignments();    

    # # set sequences to yield the reference genome sequence OUTSIDE of inverted segments
    # # remember that aligner has already RC'ed REVERSE alignments
    # # two read SEQs exit on same strand of source molecule (not of the genome reference)
    # if($refAlns[$read1][STRAND] != $refAlns[$read2][STRAND]){
    #     if($read1 == READ1){
    #         if($refAlns[$read1][STRAND]){ # ensure both SEQ on same strand of _source_ molecule, like PDL
    #             rc(\$refAlns[$read1][SEQ]);
    #             $refAlns[$read1][QUAL] = scalar reverse($refAlns[$read1][QUAL]);
    #         } else {
    #             rc(\$refAlns[$read2][SEQ]);
    #             $refAlns[$read2][QUAL] = scalar reverse($refAlns[$read2][QUAL]);  
    #         }        
    #     } else {
    #         if($refAlns[$read1][STRAND]){ # ensure both SEQ on same strand of _source_ molecule, like PDL
    #             rc(\$refAlns[$read2][SEQ]);
    #             $refAlns[$read2][QUAL] = scalar reverse($refAlns[$read2][QUAL]);
    #         } else {
    #             rc(\$refAlns[$read1][SEQ]);
    #             $refAlns[$read1][QUAL] = scalar reverse($refAlns[$read1][QUAL]);  
    #         }    
    #     }    
    # }
    
    # output the molecule summary
    print summarizeMolecule(
        $molId, $ampliconId, $nOverlapBases, $molCount, 0,
        $read1, $read2, $side1, $side2, $groupPos1, $groupPos2
    ), "\n";   
}
sub setUnmergedGroupPos {
    my ($aln, $actingRead) = @_;
    my $groupPosI = $actingIs{"GROUP_POS".$actingRead};
    my $clipI     = $actingIs{"CLIP".$actingRead};
    my $sideI     = $actingIs{"SIDE".$actingRead};
    if(!$$aln[REVERSE]){ # i.e., a forward read
        $$aln[$groupPosI] = $$aln[POS] - $$aln[$clipI];
        $$aln[$sideI] = "R";
    } else {
        $$aln[$groupPosI] = getEnd($$aln[POS], $$aln[CIGAR]) + $$aln[$clipI];
        $$aln[$sideI] = "L";
    } 
}

#---------------------------------------------------------------------
# steps common to merged and unmerged read pairs
#---------------------------------------------------------------------

# set flags for alignment orientation and inferred genome strand
sub setAlignmentStrands {

    # extract just the REVERSE bit from the SAM FLAG for each read
    $refAlns[READ1][REVERSE] = ($refAlns[READ1][FLAG] & _REVERSE) >> SHIFT_REVERSE;  # 0=forward, 1=reverse
    $refAlns[READ2][REVERSE] = ($refAlns[READ2][FLAG] & _REVERSE) >> SHIFT_REVERSE;
    
    # use the REVERSE bit to label the reference STRAND (may modify later)
    $refAlns[READ1][STRAND] = $refAlns[READ1][REVERSE]; 
    $refAlns[READ2][STRAND] = $refAlns[READ2][REVERSE]; # may be changed later when unmerged
}

# sort reads so top and bottom strands of source duplex molecule have matching signatures for grouping
sub sortReferenceAlignments {
    my $ra1 = $refAlns[READ1];
    my $ra2 = $refAlns[READ2];
    if($$ra1[STRAND] == $$ra2[STRAND]){ # PDL and some T
        return $$ra1[STRAND] ? @botStrandIs : @topStrandIs; 
    } elsif($$ra1[RNAME] ne $$ra2[RNAME]){ # T, sort by chrom
        return $$ra1[RNAME] gt $$ra2[RNAME] ? @botStrandIs : @topStrandIs;
    } else { # I, sort by pos
        return $$ra1[POS] > $$ra2[POS] ? @botStrandIs : @topStrandIs;
    }   
}

# molecule-defining data bits
sub summarizeMolecule {
    my ($molId, $ampliconId, $nOverlapBases, $molCount, $isMerged,
        $read1, $read2, $side1, $side2, $groupPos1, $groupPos2) = @_;  
    my $nAlns1 = scalar(@{$alns[$read1]});
    my $nAlns2 = $isMerged ? 0 : scalar(@{$alns[$read2]});
    my $nAlns  = $nAlns1 + $nAlns2;
    my $nBases = $isMerged ? 
        length(${$refAlns[$read1]}[SEQ]) :
        $nOverlapBases eq "NA" ?
            "NA" :
            length(${$refAlns[$read1]}[SEQ]) + length(${$refAlns[$read2]}[SEQ]) - $nOverlapBases;
    my $span = ${$refAlns[$read2]}[$groupPos2] - ${$refAlns[$read1]}[$groupPos1] + 1;
    join("\t", 
        $ampliconId, $molCount, $isMerged ? "TRUE" : "FALSE", $nOverlapBases,
        $nAlns1, $nAlns2, $nAlns, $isMerged ? $nAlns - 1 : $nAlns - 2,
        @{$refAlns[$read1]}[RNAME, $side1, $groupPos1], # node 1 
        @{$refAlns[$read2]}[RNAME, $side2, $groupPos2], # node 2
        ${$refAlns[$read1]}[POS],
        ${$refAlns[$read2]}[POS],
        $nBases, $span
        # join("," 
        #     map()
        # )
    );
}
