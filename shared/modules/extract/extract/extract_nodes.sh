# action:
#     extract SV nodes and molecule spans from name-sorted bam
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
#     $EXTRACT_PREFIX
# input:
#     $NAME_BAM_FILE (possibly overridden by $BAM_FILE)
# outputs:
#     $COORDINATE_BAM_FILE
#     node and span files from extract_nodes.pl

# log file feedback
echo "extracting SV data from name-sorted bam"
source $GENOMEX_MODULES_DIR/source/check_name_bam_file.sh
echo $NAME_BAM_FILE

# extract SV info from name sorted bam and send all reads to coordinate sort
slurp -s 250M $NAME_BAM_FILE |
samtools view - |
perl $ACTION_DIR/extract/extract_nodes.pl |
checkPipe
# tee >( 
#     samtools view - |
#     perl $ACTION_DIR/extract/extract_nodes.pl
# ) |
# samtools sort $CRAM_OUTPUT_OPTIONS --threads $N_CPU -m $SORT_RAM_PER_CPU_INT -T $TMP_FILE_PREFIX.samtools.sort - |
# slurp -s 250M -o $COORDINATE_BAM_FILE
# checkPipe

# TODO: index coordinate bam

# plot the insert size distribution(s)
echo "plotting insert size histogram"
zcat $EXTRACT_PREFIX.insert_sizes.*.gz |
sort -k1,1n |
bedtools groupby -g 1 -c 2,3,4,5,6 -o sum | 
Rscript $ACTION_DIR/extract/insertSizes.R
checkPipe

# clean up
rm -rf $TMP_DIR_WRK/*

echo "done"
