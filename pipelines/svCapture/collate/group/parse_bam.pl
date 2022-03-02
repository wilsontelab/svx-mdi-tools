use strict;
use warnings;

# extract information that uniquely identifies each input DNA molecule
# when possible merge additional read pairs with alignment guidance
# output is one line per read pair, suitable for subsequent sorting and grouping

# initialize reporting
our $script = "parse_bam";
my ($nInputAlns, $nReadPairs, $nQualityPairs) = (0) x 20;

# load dependencies
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general IUPAC smith_waterman);
resetCountFile();

# environment variables
fillEnvVar(\my $IS_USER_BAM,  'IS_USER_BAM');
fillEnvVar(\my $N_CPU,        'N_CPU');
fillEnvVar(\my $MIN_MAPQ,     'MIN_MAPQ');
fillEnvVar(\my $LIBRARY_TYPE, 'LIBRARY_TYPE');
fillEnvVar(\my $MIN_CLIP,     'MIN_CLIP');
fillEnvVar(\my $MIN_MERGE_OVERLAP, 'MIN_MERGE_OVERLAP');
fillEnvVar(\my $MIN_MERGE_DENSITY, 'MIN_MERGE_DENSITY');

# working variables
my $collapseStrands = ($LIBRARY_TYPE eq 'TruSeq'); # otherwise, Nextera, where same signatures on opposite strands are unique source molecules

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
    #-------------
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
    _IS_PAIRED => 1, # SAM FLAG bits
    _REVERSE => 16,
    _SECOND_IN_PAIR => 128,
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
our $maxShift = 3; # for Smith-Waterman
my $nullSymbol = '*';
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
    if($threadName and $qName ne $threadName){
        $nReadPairs++;        
        print $writeH END_READ_PAIR, "\n";
        $writeH = $writeH[$nReadPairs % $N_CPU + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $qName;
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
    my ($molId, $umi1, $umi2, $isMerged);
    while(my $line = <$readH>){
        chomp $line;
        if($line eq END_READ_PAIR){
            
            # extract information on the source read pair
            if($IS_USER_BAM){ # user-supplied bam file cannot support UMIs, get merge status from FLAG
                $umi1 = $umi2 = 1;            
                $isMerged = !($alns[READ1] ? $alns[READ1][0][FLAG] & _IS_PAIRED : 1);
            } else { # reads were aligned by genomex-mdi-tools align
                ($molId, $umi1, $umi2, $isMerged) = $alns[READ1] ? split(":", $alns[READ1][0][QNAME]) : ();
            }
 
            # discard orphan reads since cannot associate UMIs with endpoints
            # process others according to current merge state
            if($isMerged or ($alns[READ1] and $alns[READ2])){
                @umis = ($umi1, $umi2);            
                $isMerged ? processPremerged(): processUnmerged();
            }

            # prepare for next read-pair
            @alns = ();    
            
        } else { # add new alignment to growing read pair
            my @aln = (split("\t", $line, ALIGN_SPLIT))[QNAME..QUAL];
            push @{$alns[($aln[FLAG] & _SECOND_IN_PAIR) >> SHIFT_SECOND]}, \@aln; # 0=read1, 1=read2
        }
    }
    
    # print quality pairs for this thread
    printCount($nQualityPairs, "nQualityPairs_$childN", "read pairs passed MAPQ $MIN_MAPQ in thread $childN");
}

#---------------------------------------------------------------------
# handle read pairs that were pre-merged by fastp
#---------------------------------------------------------------------
sub processPremerged {
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
        
        # but reject low MAPQ outermost alignments as reference alignments                       
        my $i = 0;
        while($sorted[$i][MAPQ] < $MIN_MAPQ and $i < $#sorted){ $i++ }
        $refAlns[READ1] = $sorted[$i];        
        $i = $#sorted;
        while($sorted[$i][MAPQ] < $MIN_MAPQ and $i > 0){ $i-- }
        $refAlns[READ2] = $sorted[$i];  
    }
    
    # if no alignments were of sufficient MAPQ, reject the merged read
    $refAlns[READ1][MAPQ] >= $MIN_MAPQ or return;
    $nQualityPairs++;

    # set strandedness information
    setAlignmentStrands();

    # set hypothetical outermost mapped position as used for read grouping
    setMergedGroupPos($refAlns[READ1], READ1);
    setMergedGroupPos($refAlns[READ2], READ2);
    
    # order duplex endpoints to set molecule strand
    my ($read1, $read2, $groupPos1, $groupPos2, $clip1, $clip2, $side1, $side2) = sortReferenceAlignments();

    # determine which alignment yields the reference genome sequence OUTSIDE of inverted segments
    # remember that aligner has already RC'ed REVERSE alignments
    my $seqRead = (
        $refAlns[$read1][STRAND] == $refAlns[$read2][STRAND] or
        ($read1 == READ1 and !$refAlns[READ1][STRAND]) or
        ($read1 == READ2 and  $refAlns[READ2][STRAND])
    ) ? $read1 : $read2;

    # output all fields required for read-pair consensus calling into one sortable line
    print join("\t",              
        getMoleculeKey($read1, $read2, $clip1, $clip2, $side1, $side2),
        $refAlns[$read1][$groupPos1], $refAlns[$read2][$groupPos2],
        @{$refAlns[$seqRead]}[SEQ, QUAL], ($nullSymbol) x 2,
        $refAlns[READ1][QNAME], # QNAME used for read recovery from primary bam
        $read1, # the strand of source molecule, 0 or 1 (1 means we flipped the reads)
        2, int(rand(10)) # sortable merge status
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
            
            # but reject low MAPQ outermost alignments as reference alignments
            # low MAPQ alignments are more likely to be placed differently by aligner
            # even when they are longer (e.g. randomly choosing between genome repeat units)                        
            my $i = 0;
            while($sorted[$i][MAPQ] < $MIN_MAPQ and $i < $#sorted){ $i++ }
            push @refAlns, $sorted[$i];
        } 
    }
    
    # if either alignment was not of sufficient MAPQ, reject the read pair
    ($refAlns[READ1][MAPQ] >= $MIN_MAPQ and $refAlns[READ2][MAPQ] >= $MIN_MAPQ) or return;
    $nQualityPairs++;
    
    # set strandedness information; REVERSE as aligned, STRAND corrected to source molecule
    setAlignmentStrands();
    $refAlns[READ2][STRAND] = ($refAlns[READ2][STRAND] + 1) % 2;

    # set hypothetical outermost mapped position as used for read grouping
    setUnmergedGroupPos($refAlns[READ1], READ1);
    setUnmergedGroupPos($refAlns[READ2], READ2);
    
    # set the position offset used for subsequent alignment-assisted merging
    my $end1Offset = getEnd1Offset($end1Alns[READ1], $end1Alns[READ2]);

    # merge if possible otherwise commit non-overlapping unmerged pairs
    defined($end1Offset) ? mergeWithGuidance($end1Offset) : commitUnmerged();     
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

# get the alignment position difference between two alignments on the same side of a source molecule
# applicable to paired, i.e., unmerged reads when they have sufficient overlap to generate a proper pair at one end
sub getEnd1Offset { 
    my ($outerAln1, $innerAln2) = @_;
    if($$outerAln1[RNAME] ne $$innerAln2[RNAME]){ # translocation
        return; # not proper at same end of molecule, presumably reads do not overlap sufficiently
    } else{
        my $rev1 = $$outerAln1[FLAG] & _REVERSE; # REV not necessarily set on both alns yet
        my $rev2 = $$innerAln2[FLAG] & _REVERSE;
        if($rev1 == $rev2){ # inversion
            return; # again, presume insufficient overlap, no merge guidance
        } else {
            length($$outerAln1[SEQ]) != length($$innerAln2[SEQ]) and return; # inconsistent adapter clip
            my $clip11 = $rev1 ? $rightClip_ : $leftClip_;
            $clip11 = $$outerAln1[CIGAR] =~ m/$clip11/ ? $1 : 0;
            my $clip12 = $rev2 ? $leftClip_ : $rightClip_;
            $clip12 = $$innerAln2[CIGAR] =~ m/$clip12/ ? $1 : 0;
            my $end1Offset = $rev1 ?
                ($$outerAln1[POS] + $clip11) - ($$innerAln2[POS] + $clip12) : # reverse
                ($$innerAln2[POS] - $clip12) - ($$outerAln1[POS] - $clip11);  # forward
            if(abs($end1Offset) >= length($$outerAln1[SEQ])){
                return; # reads do not overlap
            }           
            return $end1Offset; # could be negative, zero or positive
        }
    }        
}

# alignments indicate that reads do overlap, work to merge them
# can be slow to fit each molecule, but is important to merge molecules that truly do overlap
sub mergeWithGuidance {
    my ($end1Offset) = @_;
    
    # check for sufficient potential overlap
    my $readLen = length($end1Alns[READ1][SEQ]);  # not necessarily READ_LEN if adapters were trimmed
    my $overlapLen = $readLen - abs($end1Offset); # integer > 0
    $overlapLen >= $MIN_MERGE_OVERLAP or return commitUnmerged();    
    
    # extract the potential overlap sequences and qualities
    my ($refRead, $qryRead) = $end1Alns[READ1][REVERSE] ? (READ2, READ1) : (READ1, READ2);
    my ($refSeq, $qrySeq, $refQual, $qryQual);
    if($end1Offset == 0){
        $refSeq  = $end1Alns[$refRead][SEQ];
        $qrySeq  = $end1Alns[$qryRead][SEQ];
        $refQual = $end1Alns[$refRead][QUAL];
        $qryQual = $end1Alns[$qryRead][QUAL];
    } elsif($end1Offset > 0) {
        $refSeq  = substr($end1Alns[$refRead][SEQ], -$overlapLen);
        $qrySeq  = substr($end1Alns[$qryRead][SEQ], 0, $overlapLen);
        $refQual = substr($end1Alns[$refRead][QUAL], -$overlapLen);
        $qryQual = substr($end1Alns[$qryRead][QUAL], 0, $overlapLen);
    } else {
        $refSeq  = substr($end1Alns[$refRead][SEQ], 0, $overlapLen);         
        $qrySeq  = substr($end1Alns[$qryRead][SEQ], -$overlapLen);
        $refQual = substr($end1Alns[$refRead][QUAL], 0, $overlapLen);         
        $qryQual = substr($end1Alns[$qryRead][QUAL], -$overlapLen);  
    }

    # if an exact match (usually shorter overlaps skipped by fastp) accept as is
    our ($overlapSeq, $overlapQual) = ("", "");
    if($refSeq eq $qrySeq){
        $overlapSeq  = $refSeq;    
        $overlapQual = $refQual;
        
    # otherwise attempt aggressive Smith-Waterman of read against read to make them merge
    } else {
        my ($qryOnRef, $score, $startQry, $endQry, $startRef, $endRef) = smith_waterman($qrySeq, $refSeq, 1); # fast algorithm with minimal register shifting
        $score or return commitUnmerged(); # no alignment at all
        
        # prepare to track base matching
        my @ref = split('', $refSeq);
        my @qry = split('', $qrySeq);
        our $nMatches = 0;
        sub checkBaseMatch{
            my ($refBase, $qryBase) = @_;
            if($refBase eq $qryBase){ # fails on m, I or D, succeeds on M
                $nMatches++;
                $overlapSeq  .= $refBase;
                $overlapQual .= "F"; # high _sequencer_ quality since confirmed by both reads
            } else {
                $overlapSeq  .= "N";
                $overlapQual .= "!"; 
            } 
        }
        
        # add the proximal stem of ref that query did not extend into
        my $refHeadLen = max(0, $startRef - $startQry);
        if($refHeadLen > 0){
            $overlapSeq  .= substr($refSeq,  0, $refHeadLen); 
            $overlapQual .= substr($refQual, 0, $refHeadLen);
            $nMatches += $refHeadLen;
        }

        # add the bases in ref where qry was clipped and not aligned (force match base by base)
        if($startRef > 0 and $startQry > 0){
            my $leftFailLen = min($startRef, $startQry);
            my $qryOffset = $startQry - $startRef;
            foreach my $i($refHeadLen..($refHeadLen + $leftFailLen - 1)){
                checkBaseMatch($ref[$i], $qry[$i + $qryOffset]);
            }
        }        

        # assemble and add the overlap alignment
        foreach my $i(0..$#$qryOnRef){
            checkBaseMatch($ref[$i + $startRef], $$qryOnRef[$i]);
        }

        # add the bases in qry where ref was clipped and not aligned (force match base by base)
        my $maxI = $overlapLen - 1;
        if($endRef < $maxI and $endQry < $maxI){
            my $rightFailLen = min($maxI - $endRef, $maxI - $endQry);
            my $startI = $endRef + 1;            
            my $qryOffset = $endQry - $endRef;
            foreach my $i($startI..($startI + $rightFailLen - 1)){
                checkBaseMatch($ref[$i], $qry[$i + $qryOffset]);
            }
        }        

        # add the distal stem of qry that ref did not extend into
        my $qryTailLen = $endRef - $endQry;        
        if($qryTailLen > 0){
            $overlapSeq  .= substr($qrySeq,  -$qryTailLen);
            $overlapQual .= substr($qryQual, -$qryTailLen);
            $nMatches += $qryTailLen;
        }
        
        # check for sufficient merge quality
        ($nMatches >= $MIN_MERGE_OVERLAP and
         $nMatches / $overlapLen >= $MIN_MERGE_DENSITY) or return commitUnmerged();      
    }

    # assemble the merged read SEQ and QUAL
    # RC not needed since both reads are on the same strand and aligner already RC'ed them
    my ($mergeSeq, $mergeQual);
    if($end1Offset > 0){ # read pairs not also overrunning on right will fail merging and not get here
        $mergeSeq  = substr($end1Alns[$refRead][SEQ], 0, $end1Offset)
                     .$overlapSeq.
                     substr($end1Alns[$qryRead][SEQ], -$end1Offset);
        $mergeQual = substr($end1Alns[$refRead][QUAL], 0, $end1Offset)
                     .$overlapQual.
                     substr($end1Alns[$qryRead][QUAL], -$end1Offset);
    } else {
        $mergeSeq  = $overlapSeq;
        $mergeQual = $overlapQual;
    }
    
    # set additional bits for committing
    $refAlns[READ1][QNAME] =~ s/:0$/:1/; # merge status 1 is lower priority molecule than 2=fastp merge
    my ($read1, $read2, $groupPos1, $groupPos2, $clip1, $clip2, $side1, $side2) = sortReferenceAlignments(); 
    
    # output all fields required for read-pair consensus calling into one sortable line
    print join("\t",              
        getMoleculeKey($read1, $read2, $clip1, $clip2, $side1, $side2),
        $refAlns[$read1][$groupPos1], $refAlns[$read2][$groupPos2],
        $mergeSeq, $mergeQual, ($nullSymbol) x 2,
        $refAlns[READ1][QNAME], # QNAME used for read recovery from primary bam
        $read1, # the strand of source molecule, 0 or 1 (1 means we flipped the reads)
        1, int(rand(10)) # sortable merge status
    ), "\n";  
}

# reads could not be merged even with alignment guidance
# most likely non-overlapping, i.e., large source molecules
sub commitUnmerged {

    # order duplex endpoints to set molecule strand
    my ($read1, $read2, $groupPos1, $groupPos2, $clip1, $clip2, $side1, $side2) = sortReferenceAlignments();    

    # set sequences to yield the reference genome sequence OUTSIDE of inverted segments
    # remember that aligner has already RC'ed REVERSE alignments
    # two read SEQs exit on same strand of source molecule (not of the genome reference)
    if($refAlns[$read1][STRAND] != $refAlns[$read2][STRAND]){
        if($read1 == READ1){
            if($refAlns[$read1][STRAND]){ # ensure both SEQ on same strand of _source_ molecule, like PDL
                rc(\$refAlns[$read1][SEQ]);
                $refAlns[$read1][QUAL] = scalar reverse($refAlns[$read1][QUAL]);
            } else {
                rc(\$refAlns[$read2][SEQ]);
                $refAlns[$read2][QUAL] = scalar reverse($refAlns[$read2][QUAL]);  
            }        
        } else {
            if($refAlns[$read1][STRAND]){ # ensure both SEQ on same strand of _source_ molecule, like PDL
                rc(\$refAlns[$read2][SEQ]);
                $refAlns[$read2][QUAL] = scalar reverse($refAlns[$read2][QUAL]);
            } else {
                rc(\$refAlns[$read1][SEQ]);
                $refAlns[$read1][QUAL] = scalar reverse($refAlns[$read1][QUAL]);  
            }    
        }    
    }
    
    # output all fields required for read-pair consensus calling into one sortable line
    print join("\t",              
        getMoleculeKey($read1, $read2, $clip1, $clip2, $side1, $side2),
        $refAlns[$read1][$groupPos1], $refAlns[$read2][$groupPos2],
        @{$refAlns[$read1]}[SEQ, QUAL], @{$refAlns[$read2]}[SEQ, QUAL],
        $refAlns[READ1][QNAME], # QNAME used for read recovery from primary bam
        $read1, # the strand of source molecule, 0 or 1 (1 means we flipped the reads)
        0, int(rand(10)) # sortable merge status
    ), "\n";   
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

# molecule-defining data bits except for positions, which are sorted separately
sub getMoleculeKey {
    my ($read1, $read2, $clip1, $clip2, $side1, $side2) = @_;

    # include a flag whether a refAln pos was clipped, i.e., should be stratified as an SV position
    # this prevents an SV outer clip from grouping with a proper fragment with the same groupPos
    my $isSVClip1 = $refAlns[$read1][$clip1] >= $MIN_CLIP ? 1 : 0;
    my $isSVClip2 = $refAlns[$read2][$clip2] >= $MIN_CLIP ? 1 : 0;    
    
    # for Nextera libraries (or any without Y-adapter-like source strand discrimination)
    # include inferred source strand as part of the molecule key
    # since identical endpoints with opposite read1/2 orientation are independent molecules
    my $groupingStrand = $collapseStrands ? '0' : $read1; # does NOT break molecule group if TruSeq

    # return the key, consistent across merged and unmerged read pairs
    join(",", 
        $umis[$read1], @{$refAlns[$read1]}[RNAME, $side1], $isSVClip1,
        $umis[$read2], @{$refAlns[$read2]}[RNAME, $side2], $isSVClip2,
        $groupingStrand
    );
}
