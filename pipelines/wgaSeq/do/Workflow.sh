#!/bin/bash

# search all subdirectories of INPUT_DIR for cell FASTQ files
cd $INPUT_DIR
export CELL_DIRS=`ls -1 ./*/*_R1_*.fastq.gz | perl -ne '$_=~m|\./(.+)/| and print "$1\n"'`

# set derivative genome variables
source $MODULES_DIR/genome/set_genome_vars.sh

# align reads to genome and create temporary output files
runWorkflowStep 1 align align/align.sh

# make a table with cells counts per fixed-width bin
runWorkflowStep 2 count count/count.sh

# make a table of fixed-width bin x cell, plus bin metrics
runWorkflowStep 3 pivot pivot/pivot_table.sh
