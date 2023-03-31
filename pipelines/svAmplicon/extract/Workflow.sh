#!/bin/bash

# set derivative environment variables and file paths
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# set the sort parameters
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
export SORT_RAM_PER_CPU_INT=$(($RAM_PER_CPU_INT - 1000000000))
export MAX_SORT_RAM_INT=$(($TOTAL_RAM_INT - 4000000000))

# set the read length from unmerged read pairs
if [ "$READ_LEN" = "" ]; then
    export READ_LEN=`
        samtools view $NAME_BAM_FILE 2>/dev/null | 
        awk '$1 ~ /:2$/ { print length($10) }' | 
        head -n 1
    `
fi
if [ "$READ_LEN" = "" ]; then # a fallback in the very rare case that all read pairs were merged
    export READ_LEN=$((MAX_INSERT_SIZE / 4))
fi
echo "READ_LEN=$READ_LEN"

# group read pairs to unique molecule sequences, merge additional read pairs
runWorkflowStep 1 group group/group.sh

# re-align processed molecules to genome
runWorkflowStep 2 realign realign/realign.sh

# summarize information on final amplicon alignments
runWorkflowStep 3 extract extract/extract.sh

# clean up
rm -fr $TMP_DIR_WRK
# rm -f $DATA_FILE_PREFIX.interim.txt.gz
# rm -f $DATA_FILE_PREFIX.*.interim.txt.gz
# rm -f $DATA_FILE_PREFIX.discovery.txt
# rm -f $DATA_FILE_PREFIX.allowed.txt
# rm -f $DATA_FILE_PREFIX.extract.fastq.gz
