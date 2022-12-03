#!/bin/bash

echo "aligning reads to genome by cell"

# align one cell at a time, with parallelization on each cell
perl $ACTION_DIR/align/align.pl
checkPipe

# merge the final output; one file for all cells, sorted
echo "merging to single final bam file"
samtools merge -b $LOG_FILE_PREFIX.align.bam.list --threads $N_CPU - > $COORDINATE_BAM_FILE
checkPipe

# index the final bam file (to support later parallization by chrom)
echo "indexing the final bam file"
samtools index -@ $N_CPU $COORDINATE_BAM_FILE
checkPipe

echo "done"
echo
