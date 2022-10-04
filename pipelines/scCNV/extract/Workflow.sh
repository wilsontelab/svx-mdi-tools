#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh

# # set the sort parameters (RAM set in step scripts)
# source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# # clean up
# rm -fr $TMP_DIR_WRK

# list files in a directory
runWorkflowStep 1 extract extract.sh
