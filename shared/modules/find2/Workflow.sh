#!/bin/bash

# set derivative environment variables and file paths
if [[ "$ON_TARGET" =  "" || "$ON_TARGET" =  "null" ]]; then export ON_TARGET=0; fi
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh
# source $MODULES_DIR/library/set_library_vars.sh

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# find structural variants across extracted nodes and molecules
runWorkflowStep 1 find find/find_svs.sh

# clean up
rm -fr $TMP_DIR_WRK
# rm -f  $EXTRACT_PREFIX.endpoints.*.gz
# rm -f  $EXTRACT_PREFIX.insert_sizes.*.gz
# rm -f  $EXTRACT_PREFIX.nodes.*.gz
# rm -f  $EXTRACT_PREFIX.strand_counts.*.gz
# rm -f  $COMPILE_PREFIX.junction_edges.gz
# rm -f  $COMPILE_PREFIX.matchedProper.gz
# rm -f  $COMPILE_PREFIX.nodes_by_proximity.txt*
# rm -f  $COMPILE_PREFIX.outer_clips.txt*
