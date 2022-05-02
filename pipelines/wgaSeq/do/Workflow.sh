#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_wga_paths.sh

# search INPUT_DIR for cell FASTQ files
# expect either:
#   $INPUT_DIR/<cell>/*.fastq.gz (directory mode, the default) or 
#   $INPUT_DIR/<cell>_*.fastq.gz (file mode)
cd $INPUT_DIR
export CELL_IDS=`ls -1 ./*/*_R1_*.fastq.gz | perl -ne '$_=~m|\./(.+)/| and print "$1\n"'`
export INPUT_MODE="directory"
if [ "$CELL_IDS" = "" ]; then
    export CELL_IDS=`ls -1 ./*_R1_*.fastq.gz | perl -ne '$_=~m|\./(.+?)_| and print "$1\n"'`
    export INPUT_MODE="file"
    if [ "$CELL_IDS" = "" ]; then
        echo -e "could not find cell fastq files in directory:\n    $INPUT_DIR"
        exit 1
    fi
fi

# set the sort parameters (RAM set in step scripts)
source $MODULES_DIR/utilities/shell/create_temp_dir.sh

# align reads to genome and create temporary output files
runWorkflowStep 1 align align/align.sh

# make a table with cells counts per fixed-width bin
runWorkflowStep 2 count count/count.sh

# make a table of fixed-width bin x cell, plus bin metrics
runWorkflowStep 3 pivot pivot/pivot_table.sh

# TODO: perform initial NBD assessments?

# clean up
rm -fr $TMP_DIR_WRK
