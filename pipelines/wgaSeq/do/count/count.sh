#!/bin/bash

echo "counting read yields before and after filtering and grouping"

# set sort RAM
SORT_RAM=$(($TOTAL_RAM_INT - 4000000000))"b"

# collect the DUP_RATE information per cell
samtools view -q $MIN_MAPQ $COORDINATE_BAM_FILE | # output is quality filtered to reflect usable reads
cut -f 9,12,13 |
sed -e 's/XC:i://' -e 's/CB:Z://' |
tee >(
    cut -f 1,3 |
    perl $ACTION_DIR/count/makeInsSizePivot.pl > $INS_SIZES_FILE
) | 
cut -f 2,3 | 
sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM -k2,2 |
bedtools groupby -g 2 -c 1,1 -o sum,count > $COUNTS_FILE
checkPipe

# use R to add the initial read pair counts (before filtering and grouping)
Rscript $ACTION_DIR/count/count.R
checkPipe

echo "done"
echo
