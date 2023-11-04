#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh
export EXTRACT_MODULE_DIR=$MODULES_DIR/svPore/extract

# set the size units and thresholds
export USE_CHR_M=1

# set file paths, written by extract, read by analyze
if [ "$INPUT_DIR" == "NA" ]; then INPUT_DIR=${TASK_DIR}/ubam; fi
export UBAM_DIR=$INPUT_DIR
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz
export EDGES_NO_SV_FILE=$EXTRACT_PREFIX.edges.no_sv.txt.gz
export EDGES_SV_FILE=$EXTRACT_PREFIX.edges.sv.txt.gz
export EDGES_TMP_FILE=$EXTRACT_PREFIX.edges.tmp.txt.gz

# pull SV reads and representative non-SV reads from pre-aligned long reads
export EXTRACT_STEP_DIR=$EXTRACT_MODULE_DIR/extract
runWorkflowStep 1 extract $EXTRACT_STEP_DIR/extract.sh

# perform read quality assessments and analyze novel SV paths
export EXTRACT_STEP_DIR=$EXTRACT_MODULE_DIR/analyze
runWorkflowStep 2 analyze $EXTRACT_STEP_DIR/analyze.sh

rm -f $EDGES_TMP_FILE
