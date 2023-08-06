#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh
export PARSE_MODULE_DIR=$MODULES_DIR/svDJ/parse

# create tmp dir for genome loading
source $MODULES_DIR/utilities/shell/create_shm_dir.sh

# discover primers and split reads into individual amplicons
runWorkflowStep 1 parse parse.sh

# clean up
rm -fr $SHM_DIR_WRK
