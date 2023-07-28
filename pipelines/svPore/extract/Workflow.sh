#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# set the size units and thresholds
export USE_CHR_M=1

# set file paths, written by extract, read by analyze
export USAM_FILE=$DATA_FILE_PREFIX.unaligned.sam.gz
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz # from low resolution first pass
export EDGES_NO_SV_FILE=$EXTRACT_PREFIX.edges.no_sv.txt.gz
export EDGES_SV_FILE=$EXTRACT_PREFIX.edges.sv.txt.gz
export EDGES_TMP_FILE=$EXTRACT_PREFIX.edges.tmp.txt.gz

# pull SV reads and representative non-SV reads from pre-aligned long reads
runWorkflowStep 1 extract extract.sh

# perform read quality assessments and analyze novel SV paths
runWorkflowStep 2 analyze ../analyze/analyze.sh

rm -f $EDGES_TMP_FILE
