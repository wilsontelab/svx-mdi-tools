#!/bin/bash

# set and create file paths and prefixes
export EXTRACT_PREFIX=$DATA_FILE_PREFIX.extract
export COMPILE_PREFIX=$DATA_FILE_PREFIX.compile
export FIND_PREFIX=$DATA_FILE_PREFIX.find
mkdir -p $PLOTS_DIR

# set derivative environment variables
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/library/set_library_vars.sh

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# extra SV nodes and molecule spans
runWorkflowStep 1 extract extract/extract_nodes.sh

# # create a copy number map and analyze the insert size distribution
# runWorkflowStep 2 coverage coverage/coverage_map.sh

# # further processing on the extracted nodes, e.g., to establish graph edges
# runWorkflowStep 3 compile compile/compile_nodes.sh

# # find structural variants across extracted nodes and molecules
# runWorkflowStep 4 find/find_junctions.sh

# clean up
rm -rf $TMP_DIR_WRK
