#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))
export MAX_SORT_RAM_INT=$(($TOTAL_RAM_INT - 4000000000))

# set the size units and thresholds
export USE_CHR_M=1

# set file paths, written by extract, read by analyze
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz
export EDGES_NO_SV_FILE=$EXTRACT_PREFIX.edges.no_sv.txt.gz
export EDGES_SV_FILE=$EXTRACT_PREFIX.edges.sv.txt.gz
export EDGES_TMP_FILE=$EXTRACT_PREFIX.edges.tmp.txt.gz # for adapter splitting
export SEQUENCES_FILE=$EXTRACT_PREFIX.sequences.txt

# pull SV evidence from aligned long reads
runWorkflowStep 1 extract extract/extract.sh

# analyze SV junctions and molecule paths
runWorkflowStep 2 analyze analyze/analyze.sh

# clean up
rm -fr $TMP_DIR_WRK
