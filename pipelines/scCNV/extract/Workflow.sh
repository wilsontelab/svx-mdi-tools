#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh

# set alignment variables
if [[ "$BAM_DIR" = "NA" || "$BAM_DIR" = "null" || "$BAM_DIR" = "" || "$BAM_DIR" = "0" ]]; then
    export ALIGNMENT_DIR=$TASK_DIR/cram # the default, assumes user ran "scCNV align"
else 
    export ALIGNMENT_DIR=$BAM_DIR # alternatively, use can supply an input folder populated with cell-level bam/cram files
fi
if [ ! -d "$ALIGNMENT_DIR" ]; then
    echo "directory not found: $ALIGNMENT_DIR"
    exit 1
fi

# search ALIGNMENT_DIR for name-sorted cell bam/cram files
# expects: $ALIGNMENT_DIR/[DATA_NAME.][GENOME.]<cell.>[name.]<bam|cram>
cd $ALIGNMENT_DIR
export ALIGNMENT_FILE_TYPE="cram"
export ALIGNMENT_FILES=`ls -1 $ALIGNMENT_DIR/*.$ALIGNMENT_FILE_TYPE 2>/dev/null`
if [ "$ALIGNMENT_FILES" = "" ]; then
    export ALIGNMENT_FILE_TYPE="bam"
    export ALIGNMENT_FILES=`ls -1 $ALIGNMENT_DIR/*.$ALIGNMENT_FILE_TYPE 2>/dev/null`
fi
if [ "$ALIGNMENT_FILES" = "" ]; then
    echo "no bam/cram files found in: $ALIGNMENT_DIR"
    exit 1
fi

# set the sort parameters (RAM set in step scripts)
source $MODULES_DIR/utilities/shell/create_temp_dir.sh

# convert bam alignments to coverage map and SV junction files, one cell at a time
runWorkflowStep 1 extract extract.sh

# combine all cells into single composite files for the entire sample
# maintain compatibility with 10x scCNV pipeline
runWorkflowStep 2 assemble assemble.sh

# TODO: clean up interim files
