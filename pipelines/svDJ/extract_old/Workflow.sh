#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# # set the sort parameters
# source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))
# export MAX_SORT_RAM_INT=$(($TOTAL_RAM_INT - 4000000000))

# extract and group V(D)J junctions
runWorkflowStep 1 extract extract/extract.sh

# # analyze SV junctions and molecule paths
# runWorkflowStep 2 analyze analyze/analyze.sh

# # clean up
# rm -fr $TMP_DIR_WRK
