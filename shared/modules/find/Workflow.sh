#!/bin/bash

# set and create file paths and prefixes
export EXTRACT_PREFIX=$DATA_FILE_PREFIX.extract
export COMPILE_PREFIX=$DATA_FILE_PREFIX.compile
export FIND_PREFIX=$DATA_FILE_PREFIX.find
mkdir -p $PLOTS_DIR

# if a collated find, force name-sorted source file to the re-alignment bam
if [ "$IS_COLLATED" != ""  ]; then
    if [ "$USE_CRAM" = "" ]; then SUFFIX=bam; else SUFFIX=cram; fi
    export BAM_FILE=$DATA_FILE_PREFIX.$GENOME.name.realigned.$SUFFIX
fi

# set derivative environment variables
if [[ "$ON_TARGET" =  "" || "$ON_TARGET" =  "null" ]]; then export ON_TARGET=0; fi
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh

#################################
# source $MODULES_DIR/library/set_library_vars.sh
export READ_LEN=151
export MAX_TLEN=850

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

############################
# TODO: pipeline framework workflow.sh needs to include PIPELINE_NAME now that status files are multi-pipeline

# extra SV nodes and molecule spans and characterize the insert size distribution
runWorkflowStep 1 extract extract/extract_nodes.sh

# if whole-genome, create a base-level copy number map
runWorkflowStep 2 coverage coverage/coverage_map.sh

# further processing on the extracted nodes, e.g., to establish graph edges
runWorkflowStep 3 compile compile/compile_nodes.sh

# find structural variants across extracted nodes and molecules
# runWorkflowStep 4 find/find_junctions.sh

# clean up
rm -rf $TMP_DIR_WRK
