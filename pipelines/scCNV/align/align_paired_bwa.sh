#------------------------------------------------------------------
# set the product bam/cram file; abort silently if exists and not forced
#------------------------------------------------------------------
if [ "$FORCE_ALIGNMENT" = "1" ]; then
    rm -f $NAME_BAM_FILE
fi
if [ -e $NAME_BAM_FILE ]; then
    echo "bam/cram file already exists"

#------------------------------------------------------------------
# check for input sequence read files
#------------------------------------------------------------------
elif [[ "$FASTQ_FILE1" = "" || "$FASTQ_FILE2" = "" ]]; then
    echo "missing sequence read file(s); expected two paired .fastq.gz files"
    exit 1  
else
 
#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------

# use fastp for one-pass adapter trimming, read merging and quality filtering
# large numbers of threads do not improve the speed
# dup evaluation is memory intensive and not needed here
fastp \
--in1 $FASTQ_FILE1 --in2 $FASTQ_FILE2 --stdout \
--dont_eval_duplication \
--merge --include_unmerged --correction \
--html /dev/null --json /dev/null 2>/dev/null |

# align to genome using BWA; soft-clip supplementary
bwa mem -p -Y -t $N_CPU $BWA_GENOME_FASTA - 2>/dev/null |

# convert to bam/cram and add mate information while still name sorted
samtools fixmate -@ $N_CPU -m --output-fmt CRAM --reference $GENOME_FASTA - - |
slurp -s 100M -o $NAME_BAM_FILE
checkPipe

fi
