#!/bin/bash

# if a collated find, force name-sorted source file to the re-alignment bam
if [ "$IS_COLLATED" != "" ]; then
    if [[ "$USE_CRAM" = "" || "$USE_CRAM" = "0" || "$USE_CRAM" = "null" ]]; then
        SUFFIX=bam
    else 
        SUFFIX=cram
    fi
    export BAM_FILE=$DATA_FILE_PREFIX.$GENOME.name.realigned.$SUFFIX
fi

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh
source $MODULES_DIR/library/set_library_vars.sh
mkdir -p $PLOTS_DIR

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# extract SV nodes and molecule spans and characterize the insert size distribution
runWorkflowStep 1 extract extract/extract_nodes.sh

# if whole-genome, create a base-level copy number map
runWorkflowStep 2 coverage coverage/coverage_map.sh

# further processing on the extracted nodes, e.g., to establish graph edges
runWorkflowStep 3 compile compile/compile_nodes.sh

# clean up
rm -fr $TMP_DIR_WRK
# rm -f  $EXTRACT_PREFIX.insert_sizes.*.gz
# rm -f  $EXTRACT_PREFIX.nodes.*.gz 
# rm -f  $EXTRACT_PREFIX.spans.*.gz 
# rm -f  $EXTRACT_PREFIX.endpoints.*.gz
# rm -f  $EXTRACT_PREFIX.strand_counts.*.gz
# rm -f  $COMPILE_PREFIX.junction_edges.gz
# rm -f  $COMPILE_PREFIX.matchedProper.gz
# rm -f  $COMPILE_PREFIX.nodes_by_proximity.txt*
# rm -f  $COMPILE_PREFIX.outer_clips.txt*
