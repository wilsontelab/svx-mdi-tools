#!/bin/bash

# set derivative environment variables
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/library/set_library_vars.sh

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))

# group read pairs to source molecules, merge read pairs, build consensus sequences
# saves a fasta file for re-alignment
runWorkflowStep 1 group/group.sh

# # re-align processed source molecules to genome
# # saves a final bam file with unique alignments
# runWorkflowStep 2 map/remap.sh

# clean up
rm -r $TMP_DIR_WRK
