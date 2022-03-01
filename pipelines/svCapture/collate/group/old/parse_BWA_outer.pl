use strict;
use warnings;

# extract information that uniquely identifies each input DNA molecule
# output is one line per read pair, suitable for subsequent sorting and grouping

# initialize reporting
our $script  = "parse_BWA";
our $error   = "$script error";
my ($nInputAlns, $nInputReadPairs, $nQualityPairs) = (0) x 20;

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/../_workflow/workflow.pl";
require "$scriptDir/../common/utilities.pl";
map { require $_ } glob("$scriptDir/../_sequence/*.pl");
resetCountFile();

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
    REV => 6, # temp field used for parsing, RNEXT never used by script
    CLIP => 7,
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

# operating parameters
my $minMapQ = $ENV{MIN_MAPQ} || 5;
my $alnSplit = QUAL + 2;
my $isTruSeq = ($ENV{LIBRARY_TYPE} eq 'TruSeq'); # otherwise, Nextera
my $isUmi = (defined($ENV{UMI_FILE}) and $ENV{UMI_FILE} ne "");
my $minInsertSize = $ENV{MIN_INSERT_SIZE};

# regular expressions for UMI and clip handling
my $umis_      = qr/\:(\w+)\:(\w+)$/;
my $leftClip_  = qr/^(\d+)S/;
my $rightClip_ = qr/(\d+)S$/;

# process data by read pair over multiple parallel threads
launchChildThreads(\&parse_BWA);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my ($threadName);
while(my $line = <STDIN>){
    $nInputAlns++;
    my ($qName) = split("\t", $line, 2);     
    if($threadName and $qName ne $threadName){
        print $writeH END_READ_PAIR, "\n";
        $nInputReadPairs++;
        $writeH = $writeH[$nInputReadPairs % $ENV{N_CPU} + 1];
    }    
    print $writeH $line; # commit to worker thread
    $threadName = $qName;
}
print $writeH END_READ_PAIR, "\n"; # finish last read pair
$nInputReadPairs++;
finishChildThreads();

# print summary information
printCount($nInputAlns,      'nInputAlns',      'input aligned segments over all read pairs');
printCount($nInputReadPairs, 'nInputReadPairs', 'input read pairs');

# child process to parse BWA map output
sub parse_BWA {
    my ($childN) = @_;
    
    # auto-flush output to prevent buffering and ensure proper feed to sort
    $| = 1;
    
    # working variables
    my @alns = ();
    
    # run BWA input one alignment at a time
    my $readH = $readH[$childN];
    while(my $line = <$readH>){
        chomp $line;
        if($line eq END_READ_PAIR){
            if(@alns == 2){ # i.e. alignments on both reads of the pair
    
                # select the two alignments that will generate the reference coordinates for read-pair matching
                # these are the outermost alignments that pass a minimum MAPQ requirement (for reproducibility)
                my @refAlns; # the two outermost alignments
                #my @sortedAlns; # in order outermost to innermost (considered using to set molDeltaLeft on SVs)
                foreach my $read(READ1, READ2){
                    
                    # calculate outer clip lengths on all alignments
                    # use to sort alignment and to determine grouping position
                    map {
                        my $clip = ($$_[FLAG] & _REV) ? $rightClip_ : $leftClip_;
                        $$_[CLIP] = $$_[CIGAR] =~ m/$clip/ ? $1 : 0;
                    } @{$alns[$read]};                    
                    
                    # most reads have only one alignment, use it
                    if(@{$alns[$read]} == 1){
                        push @refAlns, $alns[$read][0];
                        #push @sortedAlns, $alns[$read];
                    } else {                    
                    
                        # otherwise, sort the alignment in molecule order
                        # shortest outer clip is outermost alignment on the read
                        my @sorted = sort {$$a[CLIP] <=> $$b[CLIP]} @{$alns[$read]};
                        
                        # but reject low MAPQ outermost alignments as reference alignments
                        # low MAPQ alignments are more likely to be placed differently by aligner
                        # even when they are longer (e.g. randomly choosing between genome repeat units)                        
                        my $outerI = 0;
                        while($sorted[$outerI][MAPQ] < $minMapQ and
                              $outerI < $#sorted){
                            $outerI++; # move inward to hopefully find a high quality alignment
                        }
                        push @refAlns, $sorted[$outerI];
                        #push @sortedAlns, \@sorted;
                    } 
                }

                # enforce a minimum quality expectation for the reference alignment of each read
                # supplemental alignments with low MAPQ persist, but one segment of each read must be confidently mapped
                if($refAlns[READ1][MAPQ] >= $minMapQ and
                   $refAlns[READ2][MAPQ] >= $minMapQ){
                    $nQualityPairs++;
                    
                    # extract UMIs for the source molecule
                    my @umis = $isUmi ?
                        ($refAlns[READ1][QNAME] =~ m/$umis_/) :
                        (1, 1); # for read1 and read2, respectively                    

                    # extract just the REVERSE bit from the SAM FLAG for each read
                    $refAlns[READ1][REV] = ($refAlns[READ1][FLAG] & _REV) >> SHIFT_REV;  # 0=forward, 1=reverse
                    $refAlns[READ2][REV] = ($refAlns[READ2][FLAG] & _REV) >> SHIFT_REV;
                    
                    # sort reads so top and bottom strands of source molecule have matching signatures
                    # all SEQs in a pair exit on same strand of source molecule, i.e. exactly one has been RC'ed
                    # this read ordering is to ensure grouping of reads into molecules
                    my ($outRead1, $outRead2) = (READ1, READ2);
                    if($refAlns[READ1][REV] != $refAlns[READ2][REV]){ # PDL, T-FR and T-RF, sort to 1F/2R
                        $refAlns[READ1][REV] and ($outRead1, $outRead2) = (READ2, READ1);
                    } else { # I and some T
                        if($refAlns[READ1][RNAME] ne $refAlns[READ2][RNAME]){ # T-FF and T-RR, sort by chrom
                            $refAlns[READ1][RNAME] gt $refAlns[READ2][RNAME] and ($outRead1, $outRead2) = (READ2, READ1);
                        } elsif($refAlns[READ1][POS] > $refAlns[READ2][POS]){ # I-FF and I-RR, sort by pos
                            ($outRead1, $outRead2) = (READ2, READ1);
                        }
                        $refAlns[$outRead1][REV] ? # ensure both SEQ on same strand of _source_ molecule, like PDL
                            rc(\$refAlns[$outRead1][SEQ]) :
                            rc(\$refAlns[$outRead2][SEQ]);
                            $refAlns[$outRead1][QUAL] = scalar reverse($refAlns[$outRead1][QUAL]);
                            $refAlns[$outRead2][QUAL] = scalar reverse($refAlns[$outRead2][QUAL]);
                    }                    

                    # adjust the map position to set the genome coordinate used for read grouping
                    # groupPos reflects the hypothetical outermost position with respect to each reference alignment
                    my $groupPos1 = getOutermostPos($refAlns[$outRead1]);
                    my $groupPos2 = getOutermostPos($refAlns[$outRead2]);

                    # for Nextera libraries (or any without Y-adapter-like source strand discrimination)
                    # include inferred source strand as part of the molecule key
                    # since identical endpoints with opposite read1/2 orientation are independent molecules
                    my $nexteraStrand = $isTruSeq ? '0' : $outRead1; # does NOT break molecule group if TruSeq

                    # TODO: use leftmost positions of FR left ends
                    # to determine molDelta i.e. merge offset even for SV molecules
                    # process started above, pass in our output to makeConsensus
                    # OR leave the merging process to AFTER consensus build and remap #1?
                    
                    # TODO: if trimming, TLEN should be part of key (or CLIP LEN?)
                    # NO, that's not the issue, problem is that outer clips at a junction
                    # will often have the same groupPos as proper molecules
                    # they will group with those proper an pollute their consensus
                    # and likely also be suppressed as SV evidence when clip is lost in consensus
                    
                    # output all fields required for read-pair consensus calling into one sortable line
                    print join("\t",              
                        join(",", $umis[$outRead1], @{$refAlns[$outRead1]}[RNAME, REV], # molecule key, except groupPos
                                  $umis[$outRead2], @{$refAlns[$outRead2]}[RNAME, REV], $nexteraStrand),
                        $groupPos1, $groupPos2,
                        @{$refAlns[$outRead1]}[POS, CIGAR, SEQ, QUAL], # CIGAR used when merging proper only (not I)
                        @{$refAlns[$outRead2]}[POS, CIGAR, SEQ, QUAL, QNAME], # QNAME used for read recovery from primary bam
                        $outRead1 # designates the strand of source molecule, 0 or 1 (1 means we flipped the reads)
                    ), "\n";               
                }
            }

            # prepare for next read-pair
            @alns = ();            
            
        } else { # add new alignment to growing read pair
            my @aln = (split("\t", $line, $alnSplit))[QNAME..QUAL];
            length($aln[SEQ]) >= $minInsertSize or next;
            my $read = ($aln[FLAG] & _SECOND) >> SHIFT_SECOND; # 0=read1, 1=read2
            push @{$alns[$read]}, \@aln;
        }
    }
    
    # print quality pairs for this thread
    printCount($nQualityPairs, "nQualityPairs_$childN", "read pairs that passed MAPQ $minMapQ on each read in thread $childN");
}

# get hypothetical outermost position of a read relative to a specific alignment
# reflects where POS would have been if all outer clipped bases had also aligned with no indel
sub getOutermostPos { 
    my ($aln) = @_;
    if($$aln[REV]){
        getEnd($$aln[POS], $$aln[CIGAR]) + $$aln[CLIP];        
    } else {
        $$aln[POS] - $$aln[CLIP];
    }
}


