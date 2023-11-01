use strict;
use warnings;

# as compared to other svX extractions, svPore:
#   expects PAF format as input
#   expects single-end long reads, so gaps are irrelevant
#   doesn't process outer clips, they are largely irrelevant and untrustworthy in being adapters and/or low quality bases
#   tracks strands, not sides, to allow path tracking across multiple SV junctions in a single molecule
#   has no sense of canonical strands yet, since an SV might occur in different parts/orientations of a path across multiple molecules
#   reports insertions as a positive _insertion size_, microhomology as negative (as opposed to insertion = negative _overlap size_ elsewhere)

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
    #-------------
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I", 
    #-------------
    MATCH_OP      => "M",
    NULL_OP       => "X"
};

# working variables
use vars qw($MIN_SV_SIZE 
            @alnNodes @alnMapQs @alnCigars @alnAlnQs @alnTypes @alnSizes @alnInsSizes @alnAlns
            @nodes    @mapQs    @cigars    @alnQs    @types    @sizes    @insSizes    @outAlns);   
my $minCigarSvDigits = length($MIN_SV_SIZE);

#===================================================================================================
# top level molecule parsers, called by main thread loop
#---------------------------------------------------------------------------------------------------
sub processAlignedSegment { 
    my ($aln) = @_;

    # run the CIGAR string
    # base error rate reference: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
    #   In the PAF format, column 10 divided by column 11 gives the BLAST identity, which compounds the difference for every base of each indel.
    #   The latest minimap2 at github outputs gap-compressed identity at a new de:f tag, which counts each indel as one difference. 
    # from https://lh3.github.io/minimap2/minimap2.html
    #   de	f	Gap-compressed per-base sequence divergence
    # standard PAF fields
    #   00591fa8-8f17-4540-860d-36e02dc68a08    14547   34      14535   +       chr1    248387328       145089572       145104138       14421   14607   60      
    # minimap2 tags
    #   NM:i:186        ms:i:28211       AS:i:28116      nn:i:0  tp:A:P  cm:i:2123       s1:i:12565      s2:i:50 de:f:0.0074     rl:i:2679  
    # minimap2 cigar tag, only in high-accuracy 2nd pass     
    #   cg:Z:171M1D788M2I247M9D220M1I558M2D94M1I17M1I27M1D63M1D41M1D6M1D234M5D280M1I7M5D440M4D311M1D197M2I3M1D88M2D335M8D142M1I162M1D138M1D92M2D2M2D2...
    # Dorado basecall tags, only in 2nd pass
    #   qs:i:21 ns:i:136892      ts:i:10 ch:i:1104
    my ($cigar) = ($$aln[PAF_TAGS] =~ m/cg:Z:(\S+)/);
    my ($gapCompressedError) = ($$aln[PAF_TAGS] =~ m/de:f:(\S+)/); 
    parseSvsInCigar($aln, RSTART, REND, $cigar, \&commitAlignmentEdges, 1 - $gapCompressedError); # , $$aln[N_MATCHES] / $$aln[N_BASES]

    # maintain proper 5'-3' node order on bottom strand, since svPore tracks genome paths whereas alignment always come in left-right order
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
}
sub commitAlignmentEdges {  # add continguous "A" alignment edge to molecule chain
    my ($aln, $cigar, $gapCompressedIdentity) = @_; # often includes small indels in the alignment span
    # 00ab2291-0d2b-4b87-bd1a-df484884aade    30968   35      30127   -       chr7    160567428       35199008        35229151        29851   30227   60
    # 00ab2291-0d2b-4b87-bd1a-df484884aade    30968   30193   30727   +       chr16   96330374        79690622        79691156        450     578     1    
    push @alnNodes, (
        # TODO: if aln is split by parseSvsInCigar(), QSTART and QEND still reflect the span of the total alignment on query
        # if this behavior needs to change, queryPos must be updated in parseSvsInCigar()
        join("ZZ", getSignedNode(@$aln[RNAME_INDEX, RSTART, STRAND], 1), $$aln[$$aln[STRAND] eq "+" ? QSTART : QEND]), # flip Q since PAF has both query and ref in ascending, left-to-right order
        join("ZZ", getSignedNode(@$aln[RNAME_INDEX, REND,   STRAND], 0), $$aln[$$aln[STRAND] eq "+" ? QEND : QSTART])  # (note the reverse of bottom-strand alignment sets above)
    );  
    push @alnMapQs,    $$aln[MAPQ];
    push @alnCigars,   $cigar; # CIGAR is not reversed, i.e., matches rc of read if - strand      
    push @alnAlnQs,    $gapCompressedIdentity;       
    push @alnTypes,    ALIGNMENT; 
    push @alnSizes,    $$aln[REND] - $$aln[RSTART]; # alignments carry nRefBases in eventSize
    push @alnInsSizes, "NA"; # insSize not applicable for alignments
    push @alnAlns,     $aln;
}
#---------------------------------------------------------------------------------------------------
sub processSplitJunction { 
    my ($aln1, $aln2) = @_;
    my $refPos1 = $$aln1[STRAND] eq "+" ? $$aln1[REND] : $$aln1[RSTART] + 1;
    my $refPos2 = $$aln2[STRAND] eq "-" ? $$aln2[REND] : $$aln2[RSTART] + 1;
    my $jxnType = getJxnType($aln1, $aln2, $refPos1, $refPos2);  
    {
        jxnType => $jxnType,
        svSize  => getSvSize($jxnType, $refPos1, $refPos2),
        insSize => $$aln2[QSTART] - $$aln1[QEND],  # i.e., microhomology is a negative number for svPore
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
sub getJxnType {
    my ($aln1, $aln2, $refPos1, $refPos2) = @_;  
    $$aln1[RNAME_INDEX] != $$aln2[RNAME_INDEX] and return TRANSLOCATION;
    $$aln1[STRAND] ne $$aln2[STRAND] and return INVERSION;
    my $dist = $$aln1[STRAND] eq "+" ? 
        $refPos2 - $refPos1 : 
        $refPos1 - $refPos2;  
    $dist <= 0 and return DUPLICATION;
    $dist  > 1 and return DELETION; # pos delta is 1 for a continuous proper alignment
    return UNKNOWN; # should not occur, aligner should have made a contiguous alignment, consider it an error condition
}
sub getSvSize { # always a positive integer, zero if NA
    my ($jxnType, $refPos1, $refPos2) = @_;
    ($jxnType eq TRANSLOCATION or
     $jxnType eq UNKNOWN or 
     $jxnType eq ALIGNMENT) and return 0;
    abs($refPos2 - $refPos1);
}
#===================================================================================================

#===================================================================================================
# check whether a molecule is consistent with a duplex molecule read as a single-molecule foldback inversion
# if yes, only keep the first half up to (but not including) the inversion junction and mark the molecule as duplex by setting a foldback flag
# it is most sensitive and acceptable to reject any inversion with reverse-complement overlap between its flanking alignments
#   ----->
#         | V
#   <-----
# strictly speaking we expect symmetry, but cannot count on complete alignment of both flanks
# we might expect adapters in the junction, but they may not be there depending on signal quality and basecalling
# and foldback inversions are so suspect as artifacts that we maintain high sensivitiy for purging them
#---------------------------------------------------------------------------------------------------
# the resulting single read is derived from the proximal part of the source read, i.e., the first sequenced strand, which is often longer and higher quality
# that read has not been subjected to stereo basecalling, so persists as dx:i:0, foldback=TRUE
#---------------------------------------------------------------------------------------------------
# this process is independent of duplex molecules where the two strands gave two reads, which are processed and flagged by Dorado duplex
# also, as of 0.4.0, Dorado does stereo basecalling on foldback chimeras, but ONLY IF they can be found by adapter splitting
# see:  https://github.com/nanoporetech/dorado/issues/443
#---------------------------------------------------------------------------------------------------
sub checkForDuplex {
    my ($nEdges) = @_;
    $nEdges > 1 or return 0; # i.e., one stranded non-SV = not a foldback
    for(my $i = 0; $i <= $#types; $i++){
        $types[$i] eq INVERSION or next;
        if($outAlns[$i - 1][RSTART] <= $outAlns[$i + 1][REND] and
           $outAlns[$i + 1][RSTART] <= $outAlns[$i - 1][REND]){
            my $j = $i - 1;
            @mapQs  = @mapQs[0..$j]; # no need to splice @nodes, they'll never print
            @cigars = @cigars[0..$j];
            @alnQs  = @alnQs[0..$j];
            @types  = @types[0..$j];
            @sizes  = @sizes[0..$j];
            @insSizes = @insSizes[0..$j];
            return 1; # i.e., two strands = a foldback duplex (may still have an SV upstream of the foldback inversion)
        }
    }
    return 0; # i.e., one stranded SV = not a foldback
}
#===================================================================================================

1;
