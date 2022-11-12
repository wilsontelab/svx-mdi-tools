#!/bin/bash

# set derivative genome variables
export GENOMEX_SUITE=genomex-mdi-tools
export GENOMEX_MODULES_DIR=$SUITES_DIR/$GENOMEX_SUITE/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh

# set alignment variables
export ALIGNMENT_DIR=$TASK_DIR/cram
mkdir -p $ALIGNMENT_DIR

# search FASTQ_DIR for cell FASTQ files
# expect either:
#   $FASTQ_DIR/<cell>/*.fastq.gz (directory mode, the default) or 
#   $FASTQ_DIR/<cell>_*.fastq.gz (file mode)
cd $FASTQ_DIR
export CELL_IDS=`ls -1 ./*/*_R1_*.fastq.gz 2>/dev/null | perl -ne '$_=~m|\./(.+)/| and print "$1\n"'`
export INPUT_MODE="directory"
if [ "$CELL_IDS" = "" ]; then
    export CELL_IDS=`ls -1 ./*_R1_*.fastq.gz 2>/dev/null | perl -ne '$_=~m|\./(.+?)_| and print "$1\n"'`
    export INPUT_MODE="file"
    if [ "$CELL_IDS" = "" ]; then
        echo -e "could not find cell fastq files in directory:\n    $FASTQ_DIR"
        exit 1
    fi
fi

# discover cells and align them to genome one at a time, with parallelization
runWorkflowStep 1 align align.sh
