#!/bin/bash

# set derivative environment variables
source $MODULES_DIR/scan/set_genome_vars.sh
source $MODULES_DIR/scan/set_read_file_vars.sh
source $MODULES_DIR/scan/set_alignment_vars.sh

# align reads to genome and create temporary output files
runWorkflowStep 1 align align_paired_bwa.sh

# # make a table with cells counts per fixed-width bin
# runWorkflowStep 2 count count.sh

# # make a table of fixed-width bin x cell, plus bin metrics
# runWorkflowStep 3 pivot pivot_table.sh
