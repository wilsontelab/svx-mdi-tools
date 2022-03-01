#!/bin/bash

# set derivative environment variables
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh

# #########################
source $MODULES_DIR/library/set_library_vars.sh
# export READ_LEN=151
# export MAX_TLEN=862

# set and create file paths and prefixes
export CONSENSUS_PREFIX=$DATA_FILE_PREFIX.consensus

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))
export MAX_SORT_RAM_INT=$(($RAM_PER_CPU_INT - 4000000000))

# group read pairs to source molecules, merge read pairs, build consensus sequences
runWorkflowStep 1 group group/group.sh

# re-align processed source molecules to genome
runWorkflowStep 2 realign realign/realign.sh

# clean up
rm -r $TMP_DIR_WRK
