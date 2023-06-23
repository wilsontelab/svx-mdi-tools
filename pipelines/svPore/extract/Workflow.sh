#!/bin/bash

# set derivative environment variables and file paths
export PAF_FILE_TYPE=all_reads
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# set the size units and thresholds
export USE_CHR_M=1

# set file paths, written by extract, read by analyze
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz # from low resolution first pass
export EDGES_NO_SV_FILE=$EXTRACT_PREFIX.edges.no_sv.txt.gz
export QNAMES_FILE=$EXTRACT_PREFIX.target.qNames.txt # for re-analysis at high resoultion

# pull SV reads and representative non-SV reads from pre-aligned long reads (any accuracy, CIGAR not required)
runWorkflowStep 1 extract extract.sh
