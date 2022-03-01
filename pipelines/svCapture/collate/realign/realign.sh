# action:
#     re-align consensus dna-seq read (pairs) to genome with bwa
# optional:
#     $USE_CRAM         [default: creates .bam file]
# input:
#     fastq files written by prior collate steps (group+consensus)
# output:
#     $NAME_REALIGNED_BAM_FILE
#     if "$USE_CRAM" != "", $NAME_REALIGNED_BAM_FILE will be .cram format

echo "re-aligning consensus reads to genome"

# re-align consensus reads to genome; soft-clip supplementary
slurp -s 250M gunzip -c $CONSENSUS_PREFIX.*.fq.gz |
bwa mem -p -Y -t $N_CPU $BWA_GENOME_FASTA - 2>$BWA_REALIGN_LOG_FILE |

# examine reads in stream to collect coverage information
perl $ACTION_DIR/realign/coverage.pl |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 250M -o $NAME_REALIGNED_BAM_FILE
checkPipe

echo "done"
echo
