#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh

# # set the sort parameters (RAM set in step scripts)
# source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# # clean up
# rm -fr $TMP_DIR_WRK

# set the interim files
# TODO: add scCNV to these names
export EXTRACT_FILE=$DATA_FILE_PREFIX.extract.RData
export ANALYZE_FILE=$DATA_FILE_PREFIX.analyze.RData
export NORMALIZE_FILE=$DATA_FILE_PREFIX.normalize.RData
# export OUTPUT_FILE=$DATA_FILE_PREFIX.output.RData

# extract and reformat data from 10x Cell Ranger DNA input
runWorkflowStep 1 extract extract/extract.sh

# analyze cells one at a time to set windows and fit to GC bias
runWorkflowStep 2 analyze analyze/analyze.sh

# perform cross-cell normalization, recall single-cell CNVs
runWorkflowStep 3 normalize normalize/normalize.sh

# # cluster cells and call CNVs shared by multiple cells
# runWorkflowStep 4 aggregate aggregate/aggregate.sh

# clean up
# TODO
