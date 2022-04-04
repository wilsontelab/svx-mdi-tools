use strict;
use warnings;

# constants
use constant {
    QNAME => 0, # SAM fields
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
    REV => 11,  # additional fields to take working values
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
    _REV => 16, # SAM FLAG bits
    _SECOND => 128,
    #-------------
    SHIFT_REV => 4, # how far to shift those bits to yield binary values
    SHIFT_SECOND => 7,
    #-------------
    READ1 => 0,
    READ2 => 1,
};

# regular expressions for clip handling
my $leftClip_  = qr/^(\d+)S/;
my $rightClip_ = qr/(\d+)S$/;
my %actingIs = (CLIP0 => CLIP1, # so e.g. "CLIP".READ1 => "CLIP1"
                CLIP1 => CLIP2,
                GROUP_POS0 => GROUP_POS1,
                GROUP_POS1 => GROUP_POS2,
                SIDE0 => SIDE1,
                SIDE1 => SIDE2);

# working variables
use vars qw($MIN_MAPQ $MIN_CLIP $isTruSeq);
# our $maxShift = 3;
# my $nullSymbol = '*';
my @topStrandIs = (READ1, READ2, GROUP_POS1, GROUP_POS2, CLIP1, CLIP2, SIDE1, SIDE2);
my @botStrandIs = (READ2, READ1, GROUP_POS2, GROUP_POS1, CLIP2, CLIP1, SIDE2, SIDE1);
my (@alns, @refAlns, @end1Alns) = ();


# # child process to parse BWA map output
# sub parse_BWA {
#     my ($childN) = @_;
    
#     # auto-flush output to prevent buffering and ensure proper feed to sort
#     $| = 1;

#     # run BWA input one alignment at a time
#     my $readH = $readH[$childN];
#     while(my $line = <$readH>){
#         chomp $line;
#         if($line eq END_READ_PAIR){
            
#             # extract information on the source read pair
#             my ($readPairId, $umi1, $umi2, $isMerged) =
#                 $alns[READ1] ? split(":", $alns[READ1][0][QNAME]) : ();
                
#             # discard orphan reads since cannot associate UMIs with endpoints
#             # process others according to current merge state
#             if($isMerged or ($alns[READ1] and $alns[READ2])){
#                 @umis = ($umi1, $umi2);            
#                 $isMerged ? processPremerged(): processUnmerged();
#             }

#             # prepare for next read-pair
#             @alns = ();    
            
#         } else { # add new alignment to growing read pair
#             my @aln = (split("\t", $line, ALIGN_SPLIT))[QNAME..QUAL];
#             push @{$alns[($aln[FLAG] & _SECOND) >> SHIFT_SECOND]}, \@aln; # 0=read1, 1=read2
#         }
#     }
    
#     # print quality pairs for this thread
#     printCount($nQualityPairs, "nQualityPairs_$childN", "read pairs passed MAPQ $minMapQ in thread $childN");
# }


#---------------------------------------------------------------------
# handle read pairs that were pre-merged by fastp
#---------------------------------------------------------------------
sub processPremerged {
    @refAlns = ();
    
    # calculate clip lengths for all alignments
    map {
        my $clip1  = ($$_[FLAG] & _REV) ? $rightClip_ : $leftClip_;
        my $clip2  = ($$_[FLAG] & _REV) ? $leftClip_ : $rightClip_;
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
        # TODO: needs chrom for coverage map   
        getMoleculeKey($read1, $read2, $clip1, $clip2, $side1, $side2),
        $refAlns[$read1][$groupPos1], $refAlns[$read2][$groupPos2],
        # @{$refAlns[$seqRead]}[SEQ, QUAL], ($nullSymbol) x 2,
        # $refAlns[READ1][QNAME], # QNAME used for read recovery from primary bam
        $read1, # the strand of source molecule, 0 or 1 (1 means we flipped the reads)
        # 2, int(rand(10)) # sortable merge status
    ), "\n";
}
sub setMergedGroupPos {
    my ($aln, $actingRead) = @_;
    my $groupPosI = $actingIs{"GROUP_POS".$actingRead};
    my $clipI     = $actingIs{"CLIP".$actingRead};
    my $sideI     = $actingIs{"SIDE".$actingRead};
    if(($actingRead == READ1 and !$$aln[REV]) or
       ($actingRead == READ2 and  $$aln[REV])){
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
            my $clip = ($$_[FLAG] & _REV) ? $rightClip_ : $leftClip_;
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
    
    # set strandedness information; REV as aligned, STRAND corrected to source molecule
    setAlignmentStrands();
    $refAlns[READ2][STRAND] = ($refAlns[READ2][STRAND] + 1) % 2;

    # set hypothetical outermost mapped position as used for read grouping
    setUnmergedGroupPos($refAlns[READ1], READ1);
    setUnmergedGroupPos($refAlns[READ2], READ2);
    
    # set the position offset used for subsequent alignment-assisted merging
    my $end1Offset = getEnd1Offset($end1Alns[READ1], $end1Alns[READ2]);

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
        # TODO: needs chrom for coverage map     
        getMoleculeKey($read1, $read2, $clip1, $clip2, $side1, $side2),
        $refAlns[$read1][$groupPos1], $refAlns[$read2][$groupPos2],
        # @{$refAlns[$read1]}[SEQ, QUAL], @{$refAlns[$read2]}[SEQ, QUAL],
        # $refAlns[READ1][QNAME], # QNAME used for read recovery from primary bam
        $read1, # the strand of source molecule, 0 or 1 (1 means we flipped the reads)
        # 0, int(rand(10)) # sortable merge status
    ), "\n";  
}
sub setUnmergedGroupPos {
    my ($aln, $actingRead) = @_;
    my $groupPosI = $actingIs{"GROUP_POS".$actingRead};
    my $clipI     = $actingIs{"CLIP".$actingRead};
    my $sideI     = $actingIs{"SIDE".$actingRead};
    if(!$$aln[REV]){ # i.e. a forward read
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
        my $rev1 = $$outerAln1[FLAG] & _REV; # REV not necessarily set on both alns yet
        my $rev2 = $$innerAln2[FLAG] & _REV;
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

#---------------------------------------------------------------------
# steps common to merged and unmerged read pairs
#---------------------------------------------------------------------

# set flags for alignment orientation and inferred genome strand
sub setAlignmentStrands {

    # extract just the REVERSE bit from the SAM FLAG for each read
    $refAlns[READ1][REV] = ($refAlns[READ1][FLAG] & _REV) >> SHIFT_REV;  # 0=forward, 1=reverse
    $refAlns[READ2][REV] = ($refAlns[READ2][FLAG] & _REV) >> SHIFT_REV;
    
    # use the REVERSE bit to label the reference STRAND (may modify later)
    $refAlns[READ1][STRAND] = $refAlns[READ1][REV]; 
    $refAlns[READ2][STRAND] = $refAlns[READ2][REV]; # may be changed later when unmerged
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

    # include a flag whether a refAln pos was clipped, i.e. should be stratified as an SV position
    # this prevents an SV outer clip from grouping with a proper fragment with the same groupPos
    my $isSVClip1 = $refAlns[$read1][$clip1] >= $MIN_CLIP ? 1 : 0;
    my $isSVClip2 = $refAlns[$read2][$clip2] >= $MIN_CLIP ? 1 : 0;    
    
    # for Nextera libraries (or any without Y-adapter-like source strand discrimination)
    # include inferred source strand as part of the molecule key
    # since identical endpoints with opposite read1/2 orientation are independent molecules
    my $nexteraStrand = $isTruSeq ? '0' : $read1; # does NOT break molecule group if TruSeq

    # return the key, consistent across merged and unmerged read pairs
    join(",", 
        @{$refAlns[$read1]}[RNAME, $side1], $isSVClip1,
        @{$refAlns[$read2]}[RNAME, $side2], $isSVClip2,
        $nexteraStrand
    );
}

1;
