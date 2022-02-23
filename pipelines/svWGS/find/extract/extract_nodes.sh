# action:
#     extract SV nodes and molecule spans from name-sorted bam
# expects:
#     source $MODULES_DIR/scan/set_genome_vars.sh
#     source $MODULES_DIR/scan/set_alignment_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
#     $EXTRACT_PREFIX
# input:
#     $NAME_BAM_FILE
# outputs:
#     $COORDINATE_BAM_FILE
#     node and span files from extract_nodes.pl

# extract SV info from name sorted bam and send all reads to coordinate sort

echo "extracting SV data from name-sorted bam"
echo $NAME_BAM_FILE
slurp -s 500M $NAME_BAM_FILE |
tee >( 
    samtools view - |
    perl $ACTION_DIR/extract/extract_nodes.pl
) |
samtools sort $CRAM_OUTPUT_OPTIONS --threads $N_CPU -m $SORT_RAM_PER_CPU_INT -T $TMP_FILE_PREFIX.samtools.sort - |
slurp -s 500M -o $COORDINATE_BAM_FILE
checkPipe

# clean up
rm -rf $TMP_DIR_WRK/*

echo "done"
echo
