# action:
#     re-align read (pairs) to genome with bwa
# optional:
#     $USE_CRAM         [default: creates .bam file]
# input:
#     fastq file written by prior collate steps
# output:
#     $NAME_REALIGNED_BAM_FILE
#     if "$USE_CRAM" != "", $NAME_REALIGNED_BAM_FILE will be .cram format

echo "re-aligning unique amplicon molecules to genome"

# re-align to genome; soft-clip supplementary
slurp -s 250M gunzip -c $COLLATE_PREFIX.fastq.gz |
bwa mem -Y -t $N_CPU $BWA_GENOME_FASTA - 2>$BWA_REALIGN_LOG_FILE | # NOT smart pairing, reads aligned individually

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m $CRAM_OUTPUT_OPTIONS - - |
slurp -s 250M -o $NAME_REALIGNED_BAM_FILE
checkPipe

echo "done"
