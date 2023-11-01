use strict;
use warnings;

# characterize the nodes of all alignments and junctions in an edge stream
# adds information used in SV filtering and analysis
# creates a coverage map with proper handling of duplex coverage, i.e., not over-counting duplexes

# load dependencies
our $script = "extend_edges";
our $error  = "$script error";
our ($matchScore, $mismatchPenalty, $gapOpenPenalty, $gapExtensionPenalty, $maxShift) = 
    (1,           -1.5,             -2.5,            -1,                   3);
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman);
resetCountFile();

# constants
use constant {
    END_READ_EDGES => "__END_READ_EDGES__",
    #=============
    QNAME_ => 0, # edge format fields
    NODE1_ => 1,
    QSTART_ => 2,    
    NODE2_ => 3,
    QEND_ => 4,
    MAPQ_ => 5,
    CIGAR_ => 6,
    GAP_COMPRESSED_IDENTITY_ => 7,
    EDGE_TYPE_ => 8,
    EVENT_SIZE_ => 9,
    INSERT_SIZE_ => 10,
    FOLDBACK_ => 11,
    #-------------
    jxnSeq => 12, # in query orientation, not reverse-complemented yet
    baseQual => 13,   
    alnBaseQual => 14, 
    alnSize => 15, # added to edges by processRead_
    # sStart => 15,
    # sEnd => 16,
    # #-------------
    # clip5 => 0, # adapter scores added to edges here by addAdaptersScores
    # score5 => 0,
    # nBases5 => 0,
    # start5 => 0,
    # end5 => 0,
    # clip3 => 0,
    # score3 => 0,
    # nBases3 => 0,
    # start3 => 0,
    # end3 => 0,
    # score5C => 0,
    # nBases5C => 0,
    # start5C => 0,
    # end5C => 0,
    # score3C => 0,
    # nBases3C => 0,
    # start3C => 0,
    # end3C => 0,
    # #-------------
    # channel => 0, # added to edges by finishEdge
    # pod5File => 0,
    # blockN => 0,
    # edgeN => 0,
    # duplex => 0,
    #=============
    _QNAME => 0, # SAM fields
    _FLAG => 1,
    _RNAME => 2,
    _POS => 3,
    _MAPQ => 4,
    _CIGAR => 5,
    _RNEXT => 6,
    _PNEXT => 7,
    _TLEN => 8,
    _SEQ => 9,
    _QUAL => 10,
    _TAGS => 11,
    #=============
    ALIGNMENT     => "A", # the single type for a contiguous aligned segment
    TRANSLOCATION => "T", # edge/junction types (might be several per source molecule)
    INVERSION     => "V",
    DUPLICATION   => "U",
    DELETION      => "D",
    UNKNOWN       => "?",
    INSERTION     => "I"
};
my %inlineTypes = map { $_ => 1 } (ALIGNMENT, DELETION, DUPLICATION, INSERTION);

# environment variables
fillEnvVar(\our $N_CPU,             'N_CPU');
fillEnvVar(\our $EXTRACT_PREFIX,    'EXTRACT_PREFIX');
fillEnvVar(\our $PIPELINE_DIR,      'PIPELINE_DIR');
fillEnvVar(\our $EXTRACT_STEP_DIR,  'EXTRACT_STEP_DIR');
fillEnvVar(\our $WINDOW_SIZE,       'WINDOW_SIZE');
fillEnvVar(\our $GENOME_FASTA,      'GENOME_FASTA');
fillEnvVar(\our $EDGES_NO_SV_FILE,  'EDGES_NO_SV_FILE');
fillEnvVar(\our $UBAM_DIR,          'UBAM_DIR');
fillEnvVar(\our $MIN_SV_SIZE,       'MIN_SV_SIZE');
fillEnvVar(\our $MIN_MAPQ,          'MIN_MAPQ');
fillEnvVar(\our $MIN_ALIGNMENT_SIZE,        'MIN_ALIGNMENT_SIZE');
fillEnvVar(\our $MIN_ALIGNMENT_IDENTITY,    'MIN_ALIGNMENT_IDENTITY');
fillEnvVar(\our $ONT_LIBRARY_TYPE,  'ONT_LIBRARY_TYPE');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# load additional dependencies
map { require "$EXTRACT_STEP_DIR/$_.pl" } qw(initialize_windows);
initializeWindowCoverage();

# open the unaligned read files for parallel threading of QNAMES
open my $ubamH, "-|", "samtools cat $UBAM_DIR/*.unaligned.bam | samtools view -" or die "could not open bam files in $UBAM_DIR: $!\n";

# open output handles
open my $nosvH,  "|-", "pigz -p $N_CPU -c | slurp -s 10M -o $EDGES_NO_SV_FILE" or die "could not open: $EDGES_NO_SV_FILE\n";

# ONT adapter information, validated here: https://github.com/rrwick/Porechop/blob/master/porechop/adapters.py
my $ADAPTER_CORE = $ONT_LIBRARY_TYPE; # if not using a standard adapters, user must provide the apapter sequence, where 3' most base is fused to gDNA at the 5' start of the read
$ONT_LIBRARY_TYPE eq "ligation" and $ADAPTER_CORE = "ACTTCGTTCAGTTACGTATTGCT";    # ligation kit, duplex portion of the adapter; last T matches the one-base A-tail
$ONT_LIBRARY_TYPE eq "rapid"    and $ADAPTER_CORE = "TTCGCGTTTTTCGTGCGCCGCTTCA";  # rapid kit, 25bp empirically determined using extract_adapters.pl
my $adapterPadding = 10;
my $ADAPTER_CORE_RC = $ADAPTER_CORE; # ADAPTER_CORE is fused to 5' genomic ends, ADAPTER_CORE_RC is fused to 3' ends
rc(\$ADAPTER_CORE_RC);
my $coreLen = length($ADAPTER_CORE);
my $outsideLen = $coreLen + $adapterPadding;
my $maxCheckLen = $coreLen + 2 * $adapterPadding;

# adapter SVM shared support functions
sub getCandidate5 {
    my ($edge, $pos5) = @_;
    $$edge[$pos5] >= $outsideLen ? {
        startPos    => $$edge[$pos5] - $outsideLen,
        length      => $maxCheckLen,
        outsideLen  => $outsideLen
     } : {
        startPos => 0,
        length => $$edge[$pos5] + $adapterPadding,
        outsideLen => $$edge[$pos5]
     }
}
sub getCandidate3 {
    my ($edge, $pos3, $readLen) = @_;
    $ONT_LIBRARY_TYPE eq "rapid" and return { # Tn5 3' ends aren't bonded to an adapter
        startPos => -1,
        length   => -1
    };
    my $endPos = $$edge[$pos3] + $outsideLen;
    my $overrun = $endPos - $readLen;
    {
        startPos => $$edge[$pos3] - $adapterPadding,
        length   => $overrun <= 0 ? $maxCheckLen : $maxCheckLen - $overrun
    }
}
sub runAdapterSW {
    my ($x, $qSeq, $ref) = @_;
    if($$x{startPos} == -1){
        $$x{sw} = {
            bestScore => 0,
            nBases    => 0,
            qryStart  => 0,
            qryEnd    => 0
        };
        return;
    }
    my $qry = substr($qSeq, $$x{startPos}, $$x{length});
    my ($qryOnRef, $score, $startQry, $endQry) = smith_waterman($qry, $ref, undef, undef, 1);
    $$x{sw} = {
        bestScore => $score,
        nBases    => $score ? scalar(@$qryOnRef) : 0,
        qryStart  => $score ? $startQry : 0,
        qryEnd    => $score ? $endQry : 0
    }
}
my @nullControl = (0) x 8;
my @nullAdapterScores = (0) x 18;
sub addAdaptersScores {
    my ($read, $edge, $isSvRead, $isJunction) = @_;
    my ($pos5, $pos3) = $isJunction ? (QEND_, QSTART_) : (QSTART_, QEND_); # edge query positions where 5' and 3' adapters are expected
    my $readLen = length($$read[_SEQ]);

    # 5' side/start of edge, expected match to adapter
    my $x5 = getCandidate5($edge, $pos5);
    runAdapterSW($x5, $$read[_SEQ], $ADAPTER_CORE);

    # 3' side/start of edge, expected match to rc(adapter)
    my $x3 = getCandidate3($edge, $pos3, $readLen);
    runAdapterSW($x3, $$read[_SEQ], $ADAPTER_CORE_RC);

    # commit the scores of the actual clip sequence of every query node
    push @$edge, (
        $$edge[$pos5],
        $$x5{sw}{bestScore}, 
        $$x5{sw}{nBases}, 
        $$x5{sw}{bestScore} > 0  ? $$x5{sw}{qryStart} - $$x5{outsideLen} : 0,
        $$x5{sw}{bestScore} > 0  ? $$x5{sw}{qryEnd}   - $$x5{outsideLen} : 0,

        $readLen - $$edge[$pos3], 
        $$x3{sw}{bestScore}, 
        $$x3{sw}{nBases},
        $$x3{sw}{bestScore} > 0 ? $$x3{sw}{qryStart} - $adapterPadding : 0,
        $$x3{sw}{bestScore} > 0 ? $$x3{sw}{qryEnd}   - $adapterPadding : 0
    );

    # for training alignments, also process paired internal control sequences
    if($isSvRead){
        push @$edge, @nullControl;
    } else {
        my $x5C = {
            startPos => $$edge[$pos5] + $maxCheckLen,
            length => $maxCheckLen,
            outsideLen => $outsideLen
        };
        runAdapterSW($x5C, $$read[_SEQ], $ADAPTER_CORE);
        my $x3C = $ONT_LIBRARY_TYPE eq "rapid" ? { # Tn5 3' ends aren't bonded to an adapter
            startPos => -1,
            length   => -1
        } : {
            startPos => $$edge[$pos3] - 2 * $maxCheckLen,
            length => $maxCheckLen
        };
        runAdapterSW($x3C, $$read[_SEQ], $ADAPTER_CORE_RC);
        push @$edge, (
            $$x5C{sw}{bestScore}, 
            $$x5C{sw}{nBases},
            $$x5C{sw}{bestScore} > 0  ? $$x5C{sw}{qryStart}  - $$x5C{outsideLen} : 0,
            $$x5C{sw}{bestScore} > 0  ? $$x5C{sw}{qryEnd}    - $$x5C{outsideLen} : 0,

            $$x3C{sw}{bestScore}, 
            $$x3C{sw}{nBases},
            $$x3C{sw}{bestScore} > 0 ? $$x3C{sw}{qryStart}  - $adapterPadding : 0,
            $$x3C{sw}{bestScore} > 0 ? $$x3C{sw}{qryEnd}    - $adapterPadding : 0
        );
    }
}

# collect all edges of a single read
# process data by read over multiple parallel threads
launchChildThreads(\&processRead);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
my $readH  = $readH[1];
my ($nReads, $nSv, $nNoSv, $nTrainable, $prevQName, @edges, @edgeArrays) = (0, 0, 0, 0);
sub printReadEdges {
    $nReads++;
    my $read = <$ubamH>;
    while($read !~ m/^$prevQName\s/){ # handle reads found in SAM that failed to align and are absent from PAF; space is important for duplex names matching!
        $read = <$ubamH>;
        $read or last;
    }
    !$read and die "$error:\ngot to end of reads in:\n$UBAM_DIR\nwithout encountering read ID:\n$prevQName\n";
    my @read = split("\t", $read, 12);
    my ($duplex)  = ($read[_TAGS] =~ m/dx:i:(\S+)/);
    my ($isSplit) =  $read[_TAGS] =~ m/pi:Z:\S+/ ? 1 : 0;
    my $isSv = (@edges > 1);
    my $isSimplex = ($duplex == 0);
    $isSv ? $nSv++ : $nNoSv++;
    my $isTrainable = (!$isSv and !$isSplit and $isSimplex and length($read[_SEQ]) > 500);
    $isTrainable and $nTrainable++;
    if($isSv or ($isTrainable and $nTrainable <= 10000)){
        foreach my $edge(@edges){
            print $writeH $edge;
            $writeH->flush();
            $readH->flush();
        }
        # print $writeH join("", @edges);
        # $writeH->flush();
        # $readH->flush();
        print $writeH END_READ_EDGES, "\t", $read;
        $writeH->flush();
        $readH->flush();
    }
    my $increment = $isSimplex ? 1 : -$duplex; # thus increment coverage up on simplex, down on duplex (thus, +1 simplex for non-duplex, +2-1 net for duplex)  
    foreach my $edgeArray(@edgeArrays){
        $$edgeArray[EDGE_TYPE_] eq ALIGNMENT or next;
        incrementWindowCoverage(@$edgeArray[NODE1_, NODE2_], $increment);
    }
    !$isSv and print $nosvH $edges[0]; # the record of all simple alignment edges
}
while(my $edge = <STDIN>){
    chomp $edge;
    my @edgeArray = split("\t", $edge);
    if($prevQName and $prevQName ne $edgeArray[QNAME_]){
        printReadEdges();
        $writeH = $writeH[$nReads % $N_CPU + 1];
        $readH  = $readH [$nReads % $N_CPU + 1];
        @edges = ();
        @edgeArrays = ();
    }    
    push @edges, $edge."\n"; # commit edge to worker thread
    push @edgeArrays, \@edgeArray;
    $prevQName = $edgeArray[QNAME_];
}
printReadEdges(); # finish last read
finishChildThreads();
printWindowCoverage();
close $ubamH;
close $nosvH;

# print summary information
printCount($nReads, 'nReads',   'total reads processed');
printCount($nSv,    'nSv',      'reads with at least one candidate SV');
printCount($nNoSv,  'nNoSv',    'single-alignment reads with no SV (up to 10K simplex kept for adapter training)');

# extend and print one molecule's edges
sub finishEdge {
    my ($edge, $duplex, $isSplit, $blockN, $edgeN, $channel, $pod5File) = @_;
    $$edge[baseQual] ne "NA" and $$edge[baseQual] = int($$edge[baseQual] * 10 + 0.5) / 10;
    push @$edge, (
        $channel || "NA",
        $pod5File || "NA",
        $blockN,
        $edgeN,
        $duplex,
        $isSplit
    );
}
sub processRead_ {
    my ($read, $edges) = @_;
    my ($duplex)   = ($$read[_TAGS] =~ m/dx:i:(\S+)/);
    my ($isSplit) =  $$read[_TAGS] =~ m/pi:Z:\S+/ ? 1 : 0;

    # process an SV, junction-containing molecule
    if(@$edges > 1){
        my ($channel)  = ($$read[_TAGS] =~ m/ch:i:(\S+)/);
        my ($pod5File) = ($$read[_TAGS] =~ m/fn:Z:(\S+)/);
        # my ($readN)  = ($$read[_TAGS] =~ m/rn:i:(\S+)/); 

        # parse qPos to sPos (i.e, sample position)
        # my @baseSamples;
        # my ($nTotalSamples) = ($$read[_TAGS] =~ m/ns:i:(\S+)/);
        # my ($nTrimmedSamples) = ($$read[_TAGS] =~ m/ts:i:(\S+)/);
        # my ($downsampling, $moves) = ($$read[_TAGS] =~ m/\tmv:B:c,(\d+),(\S+)/);
        # my @moves = split(",", $moves);
        # map { $moves[$_] and  push @baseSamples, $_ * $downsampling + $nTrimmedSamples } 0..$#moves; 

        # first run sets baseQual for each alignment and junction itself; required to set alnBaseQual below
        foreach my $i(0..$#$edges){
            my $edge = $$edges[$i];
            my $isJunction = ($i % 2);
            if(!$isJunction){
                $$edge[jxnSeq] = "NA";
                $$edge[baseQual] = getAvgQual(substr($$read[_QUAL], $$edge[QSTART_], $$edge[EVENT_SIZE_]));
            } elsif($$edge[INSERT_SIZE_] > 0){
                $$edge[jxnSeq] = substr($$read[_SEQ], $$edges[$i - 1][QEND_], $$edge[INSERT_SIZE_]);
                $$edge[baseQual] = getAvgQual(substr($$read[_QUAL], $$edges[$i - 1][QEND_], $$edge[INSERT_SIZE_]));
            } elsif($$edge[INSERT_SIZE_] < 0){
                $$edge[jxnSeq] = substr($$read[_SEQ], $$edges[$i + 1][QSTART_], -$$edge[INSERT_SIZE_]);
                $$edge[baseQual] = getAvgQual(substr($$read[_QUAL], $$edges[$i + 1][QSTART_], -$$edge[INSERT_SIZE_]));
            } else {
                $$edge[jxnSeq] = "*";
                $$edge[baseQual] = "NA";
            } 
        }
        
        # finish each edge
        my $blockN = 1;
        foreach my $i(0..$#$edges){
            my $edge = $$edges[$i];
            my $isJunction = ($i % 2);
            $$edge[alnBaseQual] = $isJunction ? min($$edges[$i - 1][baseQual],   $$edges[$i + 1][baseQual])   : "NA";
            $$edge[alnSize]     = $isJunction ? min($$edges[$i - 1][EVENT_SIZE_], $$edges[$i + 1][EVENT_SIZE_]) : "NA";
            if($isJunction and $$edge[INSERT_SIZE_] >= 5){
                addAdaptersScores($read, $edge, 1, 1); 
            } else {
                push @$edge, @nullAdapterScores;
            }

            # establish block numbering (see analyze/analyze.R for definition of blocks)
            !$inlineTypes{$$edge[EDGE_TYPE_]} and $blockN++;
            $i > 0 and !$inlineTypes{$$edges[$i - 1][EDGE_TYPE_]} and $blockN++; 
            finishEdge($edge, $duplex, $isSplit, $blockN, $i + 1, $channel, $i == 0 ? $pod5File : "NA");
        } 

    # for training molecules, just calculate adapter scores
    } else {
        my $aln = $$edges[0];
        $$aln[CIGAR_] = $$aln[jxnSeq] = $$aln[baseQual] = $$aln[alnBaseQual] = $$aln[alnSize] = "NA";
        addAdaptersScores($read, $aln); 
        finishEdge($aln, $duplex, $isSplit, 1, 1);
    }
}
sub processRead {
    my ($childN) = @_;
    my $readH = $readH[$childN];
    my @edges;
    my $tmpFile = "$EXTRACT_PREFIX.extend_edges.$childN.txt.gz";
    open my $outH, "|-", "gzip -c > $tmpFile" or die "could not open: $tmpFile: $!\n";
    while(my $line = <$readH>){
        chomp $line;
        if($line =~ m/^__END_READ_EDGES__/){
            my @read = split("\t", $line, 13);
            shift @read;
            processRead_(\@read, \@edges);
            print $outH join("\n", map { join("\t", @$_) } @edges), "\n";
            @edges = ();
        } else {
            my @edge = split("\t", $line);
            push @edges, \@edge;
        }
    }
    close $outH;
}
