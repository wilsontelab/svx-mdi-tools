# action:
#     find heterozygous and homozygous SNVs and other variants in the bulk
#     capture targets, i.e., in the source individual's constitutive genome,
#     using reads from all samples in the project, i.e., under $TASK_DIR
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $COORDINATE_BAM_FILES (one or more)
# outputs:
#     $CONSTITUTIVE_VCF, i.e., xxx.constitutive, vcf.gz

# log file feedback
echo "finding constitutive SNVs and indels in padded target regions"
echo "input bam(s): $COORDINATE_BAM_FILES"

# merge coordinate-sorted bam files into a single stream
samtools merge \
--threads $N_CPU \
-L $PADDED_TARGETS_BED \
-u -o - \
$COORDINATE_BAM_FILES |

# read pileup in target regions
bcftools mpileup \
--threads $N_CPU \
--fasta-ref $GENOME_FASTA \
--min-MQ $MIN_MAPQ \
--max-depth $MAX_READ_DEPTH \
--max-idepth $MAX_READ_DEPTH \
--ignore-RG \
--output-type u \
- |

# variant calling
bcftools call \
--threads $N_CPU \
--ploidy $PLOIDY \
--multiallelic-caller \
--variants-only \
--output-type u \
- |

# left-align and normalize indels
bcftools norm \
--threads $N_CPU \
--fasta-ref $GENOME_FASTA \
--multiallelics +any \
--output-type z \
--output $CONSTITUTIVE_VCF \
-
checkPipe

echo "done"
