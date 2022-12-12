#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh

# set output file paths
export OUTPUT_FILE=$DATA_FILE_PREFIX.$PIPELINE_NAME.normalize.rds
export PLOTS_ARCHIVE=$DATA_FILE_PREFIX.$PIPELINE_NAME.plots.tar.gz
export MANIFEST_FILE=$DATA_FILE_PREFIX.$PIPELINE_NAME.manifest.csv

# extract and reformat data from 10x Cell Ranger DNA input
# analyze cells one at a time to set windows and fit initial GC bias estimate
runWorkflowStep 1 normalize normalize.sh
