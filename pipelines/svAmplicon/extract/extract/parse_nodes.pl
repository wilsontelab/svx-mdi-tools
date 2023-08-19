use strict;
use warnings;

# as compared to other svX extractions, svAmplicon
#   has very little need for processing outer clips; those positions are fixed as the amplicon ends
#   unlike svCapture, but like svPore, tracks molecule paths, not single junctions

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
    RNAME_INDEX => 11,
    RSTART => 12, # for compatibility with PAF/svPore
    REND => 13,
    STRAND => 14,
    #-------------
    _IS_PAIRED => 1, # SAM FLAG bits
    _UNMAPPED => 4,
    _REVERSE => 16,
    _FIRST_IN_PAIR => 64,
    _SECOND_IN_PAIR => 128,
    _SUPPLEMENTARY => 2048,
    #-------------
    READ1 => 0, # for code readability
    READ2 => 1,
    MERGED_READ => 0,
    #-------------
    LEFT  => 0, # clip recovery direction for code readability
    RIGHT => 1,  
    #-------------
    LEFTWARD  => "L", # orientation of aligned source molecule relative to a mapped endpoint
    RIGHTWARD => "R",
    #-------------
    IS_PROPER => 'P', # proper and anomalous molecule classes
    IS_SV     => 'V',
    _SHIFT_REVERSE => 4, # how far to shift those bits to yield binary values
    #-------------
    _CLIP => 0, # outData and innData fields
    _POS  => 1,
    _SIDE => 2,
    _SEQ  => 3,
    _NODE => 4,
    #-------------
    _NULL => '*',   
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    MERGE_FAILURE => "M",
    PROPER        => "P"
};

# operating parameters
use vars qw($leftClip_ $rightClip_);
my @clips = ($leftClip_, $rightClip_);
my @sides = (RIGHTWARD,  LEFTWARD); # e.g., at a left end clip, the read continues righward from that point

# working variables
use vars qw($READ_LEN $MAX_INSERT_SIZE $MIN_SV_SIZE
            @alns @mol $amplicon
            @alnNodes @alnMapQs @alnCigars @alnAlnQs @alnTypes @alnSizes @alnInsSizes @alnAlns
            @nodes    @mapQs    @cigars    @alnQs    @types    @sizes    @insSizes    @outAlns);    

#===================================================================================================
# top level molecule parser, called by main thread loop, and it support subs
#---------------------------------------------------------------------------------------------------
sub parseReadAlignments {
    my ($readN, @alns) = @_;
    my @alnIs = sortReadAlignments(\@alns, LEFT, RIGHT);  
    foreach my $i(0..$#alnIs){
        my $aln = $alns[$alnIs[$i]];
        $$aln[RSTART] = $$aln[POS] - 1; # create PAF-like designations
        $$aln[REND] = getEnd($$aln[POS], $$aln[CIGAR]);
        $$aln[STRAND] = ($$aln[FLAG] & _REVERSE) ? "-" : "+";

        # add information on the split junction between two alignments, i.e., between the previous one and this one
        if($i > 0){
            my $prevAln = $alns[$alnIs[$i - 1]];
            my $jxn = processSplitJunction($prevAln, $aln);
            push @mapQs,    0;
            push @cigars,   "NA";
            push @alnQs,    0;
            push @types,    $$jxn{jxnType};            
            push @sizes,    $$jxn{svSize};
            push @insSizes, $$jxn{insSize};
            push @outAlns,  [];
        }         

        # parse this alignment    
        (@alnNodes, @alnMapQs, @alnCigars, @alnAlnQs, @alnTypes, @alnSizes, @alnInsSizes, @alnAlns) = ();
        parseSvsInCigar($aln, RSTART, REND, $$aln[CIGAR], \&commitAlignmentEdges, $readN);

        # maintain proper 5'-3' node order on bottom strand, since svPore tracks genome paths
        if($$aln[STRAND] eq "-"){
            @alnNodes    = reverse @alnNodes;
            @alnMapQs    = reverse @alnMapQs;
            @alnCigars   = reverse @alnCigars;
            @alnAlnQs    = reverse @alnAlnQs;
            @alnTypes    = reverse @alnTypes;
            @alnSizes    = reverse @alnSizes;
            @alnInsSizes = reverse @alnInsSizes;
            @alnAlns     = reverse @alnAlns; 
        }

        # commit the nodes and junctions
        push @nodes,    @alnNodes;
        push @mapQs,    @alnMapQs;
        push @cigars,   @alnCigars;
        push @alnQs,    @alnAlnQs;        
        push @types,    @alnTypes;        
        push @sizes,    @alnSizes;
        push @insSizes, @alnInsSizes;
        push @outAlns,  @alnAlns;
    }
}
sub commitAlignmentEdges{
    my ($aln, $cigar, $readN) = @_; # can include small indels in the alignment span
    # in this output of two nodes per alignment:
    #   all fields except POS are identical in node positions 1 and 2
    #   for bottom strand (STRAND == "-") alignments:
    #       POS is in reverse order, i.e., POS1 > POS2, to track the molecule path
    #       SEQ/QUAL are reverse complemented from canonical, as provided by aligner
    push @alnNodes, ( 
        join("ZZ", @$aln[RNAME, STRAND], $$aln[RSTART] + 1,  @$aln[SEQ, QUAL], $readN + 1), # svAmplicon nodes are always 1-referenced
        join("ZZ", @$aln[RNAME, STRAND, REND, SEQ, QUAL], $readN + 1)
    );
    push @alnMapQs,    $$aln[MAPQ];
    push @alnCigars,   $cigar; # CIGAR is not reversed, i.e., matches rc of read if - strand      
    push @alnAlnQs,    0;       
    push @alnTypes,    ALIGNMENT; 
    push @alnSizes,    $$aln[REND] - $$aln[RSTART]; # alignments carry nRefBases in eventSize
    push @alnInsSizes, "NA\tNA"; # insSize and appended jxnBases not applicable for alignments
    push @alnAlns,     $aln;
}
sub processSplitJunction {
    my (@alns) = @_;
    my (@innData);
    setEndpointData(\@{$innData[READ1]}, $alns[READ1], RIGHT, LEFT);
    setEndpointData(\@{$innData[READ2]}, $alns[READ2], LEFT,  RIGHT);
    my $jxnType = getJxnType($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2]);  
    my $svSize = getSvSize($innData[READ1], $innData[READ2], $jxnType);
    my ($insSize, $jxnBases) = getJxnStructure($alns[READ1], $alns[READ2], $innData[READ1], $innData[READ2]); 
    {
        jxnType  => $jxnType,
        svSize   => $svSize,
        insSize  => "$insSize\t$jxnBases" # i.e., microhomology is a negative number for svAmplicon
    }
}
#===================================================================================================

#===================================================================================================
# utility functions for characterizing SV junctions
#---------------------------------------------------------------------------------------------------
# for these purposes, consider the genome as one contiguous reference unit
#   i.e., with chromosomes conjoined in numerical index order
#       ---chr1---|---chr2---|--- etc.
#   such that ~half of translocations are handled as del/dup, half as inv
#---------------------------------------------------------------------------------------------------
# check to confirm that alignments truly do flank a SV and report its type
sub getJxnType {
    my ($aln1, $aln2, $innData1, $innData2) = @_;  

    # simplest case, translocation = alignment to different chromosomes
    # PDL vs. inversion sub-handling determined by isCanonicalStrand
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;

    # inversions determined by unexpected alignment orientiations
    $$innData1[_SIDE] eq $$innData2[_SIDE] and return INVERSION;

    # PDL are distinguish based on their position along the conjoined chromosomes
    my $dist = $$innData1[_SIDE] eq LEFTWARD ? 
                $$innData2[_POS] - $$innData1[_POS] : # LR pair, canonical PDL strand
                $$innData1[_POS] - $$innData2[_POS];  # RL pair, non-canonical strand (will be flipped later)
    $dist <= 0 and return DUPLICATION;
    $dist  > 1 and return DELETION; # pos delta is 1 for a continuous proper alignment
    return PROPER;
}
sub getSvSize { # always a positive integer, zero if NA
    my ($innData1, $innData2, $jxnType) = @_;
    ($jxnType eq TRANSLOCATION or
     $jxnType eq UNKNOWN or 
     $jxnType eq INSERTION or 
     $jxnType eq PROPER) and return 0;
    $jxnType eq INVERSION and return abs($$innData2[_POS] - $$innData1[_POS]);
    $jxnType eq DELETION and return $$innData2[_POS] - $$innData1[_POS] - 1;
    ($jxnType eq DUPLICATION or
     $jxnType eq MERGE_FAILURE) and return $$innData1[_POS] - $$innData2[_POS] + 1;
    return 0;
}
sub getJxnStructure {
    my ($aln1, $aln2, $innData1, $innData2) = @_;
    my $readPos1 = length($$aln1[SEQ]) - $$innData1[_CLIP];
    my $readPos2 = $$innData2[_CLIP] + 1;
    my $overlap = $readPos1 - $readPos2 + 1;
    my $jxnBases = "*"; # a blunt joint
    if($overlap > 0){  # microhomology
        $jxnBases = substr($$aln1[SEQ], $readPos2,      $overlap); # TODO: any value to reverse complementing inversions?
    } elsif($overlap < 0){ # inserted bases
        $jxnBases = substr($$aln1[SEQ], $readPos1 + 1, -$overlap);
    }
    (-$overlap, $jxnBases); # i.e., return insertion size
}
#===================================================================================================

#===================================================================================================
# utility functions for processing individual _alignment_ segments
#---------------------------------------------------------------------------------------------------
# order a set of alignments for a (merged) read in source molecule order
# shortest outer clip is outermost alignment on the read
sub sortReadAlignments {
    my ($alns, $fwdSide, $revSide) = @_;
    @$alns == 1 and return 0;
    my @outDatas;
    foreach my $i(0..$#$alns){
        setEndpointData(\@{$outDatas[$i]}, $$alns[$i], $fwdSide, $revSide);
    }
    return sort { $outDatas[$a][_CLIP] <=> $outDatas[$b][_CLIP] } 0..$#$alns;  
}
# collect information on one endpoint (outer or inner) of an alignment
sub setEndpointData {
    my ($data, $aln, $fwdSide, $revSide) = @_;
    my ($clip, $getEnd, $side) = ($$aln[FLAG] & _REVERSE) ?
        ($clips[$revSide], $revSide, $sides[$revSide]) :
        ($clips[$fwdSide], $fwdSide, $sides[$fwdSide]);
    $$data[_CLIP] = $$aln[CIGAR] =~ m/$clip/ ? $1 : 0; # the number of clipped bases
    $$data[_POS]  = $getEnd ? getEnd($$aln[POS], $$aln[CIGAR]) : $$aln[POS]; # the coordinate of the last aligned base
    $$data[_SIDE] = $side; # the direction the aligned read moves away from that coordinate
    $$data[_SEQ] = $$data[_CLIP] > 0 ? # the sequence of the clipped bases
        ($getEnd ?
            substr($$aln[SEQ], -$$data[_CLIP]) :
            substr($$aln[SEQ], 0, $$data[_CLIP])):
        _NULL;
    $$data[_NODE] = join("/", $$aln[RNAME], @$data[_SIDE, _POS]); # complete signature of the SV breakpoint position
}
#===================================================================================================

1;
