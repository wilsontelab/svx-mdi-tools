use strict;
use warnings;

# characterize the nodes of all alignments and junctions in an edge stream
# adds information used in SV filtering and analysis

# load dependencies
our $script = "extend_edges";
our $error  = "$script error";
our ($matchScore, $mismatchPenalty, $gapOpenPenalty, $gapExtensionPenalty, $maxShift) = 
    (1,           -1.5,             -2.5,            -1,                   3);
my $perlUtilDir = "$ENV{GENOMEX_MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman);
resetCountFile();

# constants
use constant {
    END_READ_EDGES => "__END_READ_EDGES__",
    #=============
    QNAME => 0, # edge format fields
    NODE1 => 1,
    QSTART => 2,    
    NODE2 => 3,
    QEND => 4,
    MAPQ => 5,
    CIGAR => 6,
    GAP_COMPRESSED_IDENTITY => 7,
    EDGE_TYPE => 8,
    EVENT_SIZE => 9,
    INSERT_SIZE => 10,
    N_STRANDS => 11,
    #-------------
    baseQual => 12,   
    alnBaseQual => 13, 
    alnSize => 14, # added to edges by processRead_
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
fillEnvVar(\our $EXTRACT_PREFIX,    'EXTRACT_PREFIX');
fillEnvVar(\our $USAM_FILE,         'USAM_FILE');
fillEnvVar(\our $N_CPU,             'N_CPU');
fillEnvVar(\our $MIN_SV_SIZE,       'MIN_SV_SIZE');
fillEnvVar(\our $MIN_MAPQ,          'MIN_MAPQ');
fillEnvVar(\our $MIN_ALIGNMENT_SIZE,        'MIN_ALIGNMENT_SIZE');
fillEnvVar(\our $MIN_ALIGNMENT_IDENTITY,    'MIN_ALIGNMENT_IDENTITY');

# open the usam file for parallel threading of QNAMES
open my $usamH, "-|", "zcat $USAM_FILE" or die "could not open $USAM_FILE: $!\n";

# ONT adapter information
# TODO: implement support for mosaic ends in Tn5-based rapid kit
my $ADAPTER_CORE = "ACTTCGTTCAGTTACGTATTGCT"; # duplex portion of the adapter; last T matches the one-base A-tail
my $ADAPTER_CORE_RC = $ADAPTER_CORE;          # ADAPTER_CORE is fused to 5' genomic ends, ADAPTER_CORE_RC is fused to 3' ends
rc(\$ADAPTER_CORE_RC);      
my $coreLen = length($ADAPTER_CORE);
my $adapterPadding = 10;
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
    my $endPos = $$edge[$pos3] + $outsideLen;
    my $overrun = $endPos - $readLen;
    {
        startPos => $$edge[$pos3] - $adapterPadding,
        length   => $overrun <= 0 ? $maxCheckLen : $maxCheckLen - $overrun
    }
}
sub runAdapterSW {
    my ($x, $qSeq, $ref) = @_;
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
    my ($pos5, $pos3) = $isJunction ? (QEND, QSTART) : (QSTART, QEND); # edge query positions where 5' and 3' adapters are expected
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
        my $x3C = {
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
my ($nReads, $nSv, $nNoSv, $prevQName, @edges) = (0, 0, 0);
sub printReadToThread {
    my $isSv = (@edges > 1);
    $nReads++;  
    $isSv ? $nSv++ : $nNoSv++;
    my $read = <$usamH>;
    while($read !~ m/^$prevQName/){ # handle reads found in SAM that failed to align and are absent from PAF
        $read = <$usamH>;
    }  
    $isSv or $nNoSv <= 10000 or return;
    print $writeH join("", @edges);
    print $writeH END_READ_EDGES, "\t", $read;  
}
while(my $edge = <STDIN>){
    my ($qName) = split("\t", $edge, 2);     
    if($prevQName and $qName ne $prevQName){
        printReadToThread();
        $writeH = $writeH[$nReads % $N_CPU + 1];
        @edges = ();
    }    
    push @edges, $edge; # commit edge to worker thread
    $prevQName = $qName;
}
printReadToThread(); # finish last read
finishChildThreads();
close $usamH;

# print summary information
printCount($nReads, 'nReads',   'total reads processed');
printCount($nSv,    'nSv',      'reads with at least one candidate SV');
printCount($nNoSv,  'nNoSv',    'single-alignment reads with no SV (up to 10K kept)');

# extend and print one molecule's edges
sub finishEdge {
    my ($edge, $blockN, $edgeN, $channel, $pod5File) = @_;
    $$edge[baseQual] ne "NA" and $$edge[baseQual] = int($$edge[baseQual] * 10 + 0.5) / 10;
    push @$edge, (
        $channel || "NA",
        $pod5File || "NA",
        $blockN,
        $edgeN
    );
}
sub processRead_ {
    my ($read, $edges) = @_;

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
                $$edge[baseQual] = getAvgQual(substr($$read[_QUAL], $$edge[QSTART], $$edge[EVENT_SIZE]));
            } elsif($$edge[INSERT_SIZE] > 0){
                $$edge[baseQual] = getAvgQual(substr($$read[_QUAL], $$edges[$i - 1][QEND], $$edge[INSERT_SIZE]));
            } elsif($$edge[INSERT_SIZE] < 0){
                $$edge[baseQual] = getAvgQual(substr($$read[_QUAL], $$edges[$i + 1][QSTART], -$$edge[INSERT_SIZE]));
            } else {
                $$edge[baseQual] = "NA";
            } 
        }
        
        # finish each edge
        my $blockN = 1;
        foreach my $i(0..$#$edges){
            my $edge = $$edges[$i];
            my $isJunction = ($i % 2);
            $$edge[alnBaseQual] = $isJunction ? min($$edges[$i - 1][baseQual],   $$edges[$i + 1][baseQual])   : "NA";
            $$edge[alnSize]     = $isJunction ? min($$edges[$i - 1][EVENT_SIZE], $$edges[$i + 1][EVENT_SIZE]) : "NA";
            # $$edge[sStart] = $baseSamples[$$edge[QSTART]] || "NA";
            # $$edge[sEnd]   = $baseSamples[$$edge[QEND]]   || "NA"; 
            if($isJunction and $$edge[INSERT_SIZE] >= 5){
                addAdaptersScores($read, $edge, 1, 1); 
            } else {
                push @$edge, @nullAdapterScores;
            }

            # establish block numbering (see analyze/analyze.R for definition of blocks)
            !$inlineTypes{$$edge[EDGE_TYPE]} and $blockN++;
            $i > 0 and !$inlineTypes{$$edges[$i - 1][EDGE_TYPE]} and $blockN++; 
            finishEdge($edge, $blockN, $i + 1, $channel, $i == 0 ? $pod5File : "NA");
        } 

    # for training molecules, just calculate adapter scores
    } else {
        my $aln = $$edges[0];
        $$aln[CIGAR] = "NA";
        $$aln[baseQual] = $$aln[alnBaseQual] = $$aln[alnSize] = "NA"; #  = $$aln[sStart] = $$aln[sEnd]
        addAdaptersScores($read, $aln); 
        finishEdge($aln, 1, 1);
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
 