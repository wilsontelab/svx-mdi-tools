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

# group read pairs to unique molecule sequences, merge additional read pairs
runWorkflowStep 1 group group/group.sh

# re-align processed molecules to genome
runWorkflowStep 2 realign realign/realign.sh

# clean up
rm -fr $TMP_DIR_WRK
rm -f $DATA_FILE_PREFIX.interim.txt.gz
rm -f $DATA_FILE_PREFIX.discovery.txt
rm -f $DATA_FILE_PREFIX.allowed.txt
