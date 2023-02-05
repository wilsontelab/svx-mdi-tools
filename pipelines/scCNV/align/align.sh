#!/bin/bash

# align one cell at a time, with parallelization on each cell
echo "aligning reads to genome by cell"
for CELL_ID in $CELL_IDS; do 
    echo $CELL_ID
    if [ "$INPUT_MODE" = "directory" ]; then
        export FASTQ_FILE1=./${CELL_ID}/*_R1_*.fastq.gz
        export FASTQ_FILE2=./${CELL_ID}/*_R2_*.fastq.gz
    else
        export FASTQ_FILE1=./${CELL_ID}_*_R1_*.fastq.gz
        export FASTQ_FILE2=./${CELL_ID}_*_R2_*.fastq.gz
    fi
    export DATA_GENOME_PREFIX=$ALIGNMENT_DIR/$DATA_NAME.$GENOME.$CELL_ID
    export NAME_BAM_FILE=$DATA_GENOME_PREFIX.name.cram
    source $ACTION_DIR/align_paired_bwa.sh
done

echo "done"
echo
