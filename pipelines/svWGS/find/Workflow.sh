#!/bin/bash

# set derivative environment variables
source $MODULES_DIR/scan/set_genome_vars.sh
source $MODULES_DIR/scan/set_alignment_vars.sh

# set and create file paths and prefixes
export EXTRACT_PREFIX=$DATA_FILE_PREFIX.extract
export COMPILE_PREFIX=$DATA_FILE_PREFIX.compile
export FIND_PREFIX=$DATA_FILE_PREFIX.find
mkdir -p $PLOTS_DIR

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# # discover library properties
# echo $NAME_BAM_FILE
# echo "extracting READ_LEN"
# export READ_LEN=`samtools view -f 1 $NAME_BAM_FILE | head -n 1 | awk '{print length($10)}'`
# echo READ_LEN = $READ_LEN
# echo "extracting MAX_TLEN"
# export MAX_TLEN=`samtools view $NAME_BAM_FILE | head -n 100000 | perl $ACTION_DIR/get_max_TLEN.pl`
# echo MAX_TLEN = $MAX_TLEN

export READ_LEN=151
export MAX_TLEN=862

# extra SV nodes and molecule spans
# runWorkflowStep 1 extract extract/extract_nodes.sh

# create a copy number map and analyze the insert size distribution
# runWorkflowStep 2 coverage coverage/coverage_map.sh

# further processing on the extracted nodes, e.g., to establish graph edges
runWorkflowStep 3 compile compile/compile_nodes.sh

# find structural variants across extracted nodes and molecules
runWorkflowStep 4 find/find_junctions.sh

# clean up
rm -rf $TMP_DIR_WRK
