#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $MODULES_DIR/files/set_svx_paths.sh

# group read pairs to source molecules, merge read pairs, build consensus sequences
runWorkflowStep 1 assemble assemble.sh

# clean up
