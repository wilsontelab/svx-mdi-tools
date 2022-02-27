use strict;
use warnings;

# create a consensus PHRED score for bases confirmed by multiple independent reads
# conflict bases are always assigned zero quality at any stage (PHRED 0 => !)

# see PHRED 33 scale: https://drive5.com/usearch/manual/quality_score.html

#======================================================================================

# in general at consensus prior to remap, output qualities still reflect the initial
# _sequencing_ quality, not necessarily our final confidence in the source molecule
# can later downgrade confidence in molecule sequences based on read/duplex coverage

# as an example, code will output PHRED 37 => F for singleton reads (a high confidence)
# this reflects that the base value is what the sequencer should have gotten at a cluster
# but it does not yet reflect the potential for pre-cluster errors

# it makes sense to not downgrade yet based on coverage, since it would confuse
# the BWA aligner in the remap step to see bases with low apparent _sequencing_ quality
# given that our SV pipeline maintains all read pairs regardless of coverage

# instead, at this stage we simply give a small quality _bump_ when the output
# sequence has been confirmed by multiple reads, to a max of PHRED 40

# thus, most bases after make_consesnsus will have one of three high quality values
#   37 = F = most bases in singleton read pairs (some lower if declared by sequencer)
#   38 = G = single-strand consensus base
#   40 = I = duplex consensus base
# or one null value
#   0  = ! = masked bases that failed consensus building; we reject regardless of sequencer

#======================================================================================

# use PHRED 37 for bases confirmed by merging, consistent with NovaSeq input qual
# might promote the quality of an individual base sequenced by both reads
# but never above the quality of the best input read base
sub getConsensusQual_37 {    
    my ($consensusSeq) = @_; 
    $consensusSeq or return; 
    $consensusSeq =~ tr/ACGT/F/; # F is phred 37
    $consensusSeq =~ tr/F/!/c;   # F is not an IUPAC base
    $consensusSeq;
}

# use PHRED 38 for bases confirmed by strand (not duplex) consensus
# indicates a higher degree of confidence in the sequending result
sub getConsensusQual_38 {    
    my ($consensusSeq) = @_; 
    $consensusSeq or return; 
    $consensusSeq =~ tr/ACGT/G/; # G is phred 38
    $consensusSeq =~ tr/G/!/c;   # G base yields G qual, as expected
    $consensusSeq;
}

# use phred 40 for consensus bases confirmed by both duplex reads
# thus, we promote the quality of consensus bases above the best input read quality
sub getConsensusQual_40 { 
    my ($consensusSeq) = @_;
    $consensusSeq or return;
    $consensusSeq =~ tr/ACGT/I/; # I is phred 40
    $consensusSeq =~ tr/I/!/c;   # I is not an IUPAC base
    $consensusSeq;
}

1;

