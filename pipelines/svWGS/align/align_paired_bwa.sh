# action:
#     align paired dna-seq read files to genome with read merging using fastp and bwa
# expects:
#     source $MODULES_DIR/scan/set_genome_vars.sh
#     source $MODULES_DIR/scan/set_read_file_vars.sh
#     source $MODULES_DIR/scan/set_alignment_vars.sh
#     $MIN_INSERT_SIZE
# optional:
#     $ADAPTER_SEQUENCE [default: merge-inherent trimming only]
#     $FORCE_ALIGNMENT  [default: don't overwrite NAME_BAM_FILE]
#     $USE_CRAM         [default: creates .bam file]
# input:
#     if FASTQ files are found (.fastq.gz) they are used
#     otherwise searches for SRA (.sra) files that are converted to FASTQ in a stream
# output:
#     $NAME_BAM_FILE
#     if "$USE_CRAM" != "", $NAME_BAM_FILE will be .cram format

#------------------------------------------------------------------
# set the product bam/cram file; abort silently if exists and not forced
#------------------------------------------------------------------
if [[ "$FORCE_ALIGNMENT" != "" && -e $NAME_BAM_FILE ]]; then
    echo "forcing overwrite of bam/cram file: $NAME_BAM_FILE"
    rm -f $NAME_BAM_FILE
fi
if [ -e $NAME_BAM_FILE ]; then
    echo "bam/cram file already exists"

#------------------------------------------------------------------
# check for input sequence read files
#------------------------------------------------------------------
elif [[ "$FASTQ_FILE1" != "" && "$FASTQ_FILE2" = "" ]]; then
    echo "missing fastq file(s); expected paired .fastq.gz files"
    exit 1
elif [[ "$FASTQ_FILE1" = "" && "$SRA_FILES" = "" ]]; then
    echo "missing sequence read file(s); expected two paired .fastq.gz files or a set of .sra files"
    exit 1  
else

#------------------------------------------------------------------
# initialize the alignment process
#------------------------------------------------------------------

# provide log feedback
echo "aligning paired end reads to genome $GENOME with pre-alignment merging"
echo "  read length = $READ_LEN"
if [ "$FASTQ_FILE1" = "" ]; then
    echo "  SRA files = $SRA_FILES" 
else
    echo "  read #1 = $FASTQ_FILE1"
    echo "  read #2 = $FASTQ_FILE2"    
fi

# if requested, perform "extra" adapter trimming (in addition to merge-inherent trimming)
if [[ "$ADAPTER_SEQUENCE" != "NA" && "$ADAPTER_SEQUENCE" != "null" && "$ADAPTER_SEQUENCE" != "" ]]; then
    ADAPTER_SEQUENCE="--adapter_sequence $ADAPTER_SEQUENCE"
else
    ADAPTER_SEQUENCE=""
fi
    
#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------

# pull reads from various source types to a consistent interleaved format
perl $ACTION_DIR/prepare_reads.pl |

# use fastp for one-pass adapter trimming, read merging and quality filtering  # --disable_quality_filtering
# large numbers of threads do not improve the speed; dup evaluation is memory intensive and not needed here
fastp \
--stdin --interleaved_in --stdout \
--dont_eval_duplication \
--length_required $MIN_INSERT_SIZE $ADAPTER_SEQUENCE \
--merge --include_unmerged --correction \
--html $FASTP_LOG_PREFIX.html --json $FASTP_LOG_PREFIX.json \
--report_title \"$DATA_NAME\" 2>/dev/null |

# align to genome using BWA
bwa mem -p -Y -t $N_CPU $BWA_GENOME_FASTA - 2>$BWA_LOG_FILE |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 100M -o $NAME_BAM_FILE
checkPipe

echo "done"

fi
