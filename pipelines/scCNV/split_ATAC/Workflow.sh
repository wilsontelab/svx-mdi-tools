#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh

# validate the expected input files
export CELL_RANGER_BAM_FILE=$CELL_RANGER_DIR/possorted_bam.bam
export CELL_RANGER_CELL_FILE=$CELL_RANGER_DIR/singlecell.csv
if [ ! -f $CELL_RANGER_BAM_FILE ]; then
    echo "missing file: $CELL_RANGER_BAM_FILE"
    exit 1
fi
if [ ! -f $CELL_RANGER_CELL_FILE ]; then
    echo "missing file: $CELL_RANGER_CELL_FILE"
    exit 1
fi

# set alignment variables
export ALIGNMENT_DIR=$TASK_DIR/bam
mkdir -p $ALIGNMENT_DIR

# use CellRanger ATAC inputs to parse good cells to individual bam files for extraction
runWorkflowStep 1 split_ATAC split_ATAC.sh


# https://support.10xgenomics.com/single-cell-atac/software/overview/welcome
# -rw-r--r-- 1 wilsonte wilsonte_root  80G Mar  3 20:02 possorted_bam.bam

# extract expects:
#   name-sorted
#   $ALIGNMENT_DIR/[DATA_NAME.][GENOME.]<cell.>[name.]<bam|cram>
#   any QNAME format is fine
