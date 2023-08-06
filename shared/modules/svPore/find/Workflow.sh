#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh
export FIND_MODULE_DIR=$MODULES_DIR/svPore/find

# set the size units and thresholds
export USE_CHR_M=1

# determine if this is a single-sample find or a multi-sample compare
echo "setting find mode"
export FIND_MODE=find
export EDGE_FILES=`ls -1 $ANALYZE_PREFIX.edges.rds 2>/dev/null`
if [ "$EDGE_FILES" = "" ]; then
    export FIND_MODE=compare
    export EDGE_FILES=`ls -1 $ANALYZE_GLOB_PREFIX.edges.rds 2>/dev/null`
    if [ "$EDGE_FILES" = "" ]; then
        echo "could not find analyze.edges.rds file(s) in:"
        echo "    $TASK_DIR"
        exit 1
    fi
fi
echo "  $FIND_MODE"

# pull SV reads and representative non-SV reads from pre-aligned long reads (high accuracy, CIGAR required)
runWorkflowStep 1 find $FIND_MODULE_DIR/find_svs.sh
