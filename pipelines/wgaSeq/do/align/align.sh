#!/bin/bash

echo "aligning reads to genome by cell"

# temp directory for sorting
export TMP_DIR_WRK=$TMP_DIR_LARGE/$PIPELINE_NAME"_align_"$DATA_NAME
mkdir -p $TMP_DIR_WRK

# align one cell at a time, with parallelization on each cell
perl $ACTION_DIR/align/align.pl
checkPipe

# merge the final output; one file for all cells, sorted
echo "merging to single final bam file"
samtools merge -b $LOG_FILE_PREFIX.align.bam.list --threads $N_CPU - > $BAM_FILE
checkPipe

# index the final bam file
echo "indexing the final bam file"
samtools index -@ $N_CPU $BAM_FILE
checkPipe

# clean up
rm -r $TMP_DIR_WRK

echo "done"
echo
