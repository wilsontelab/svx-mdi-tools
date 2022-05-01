#!/bin/bash

echo "filtering and counting reads to create a bin (rows) x cell (columns) pivot table"

# set files
BINS_DIR=$GENOME_BINS_DIR/fixed_width_bins
export BINS_FILE=$BINS_DIR/$GENOME.bins.size_$BIN_SIZE.k_$KMER_LENGTH.e_$N_ERRORS.bed.gz
if [ ! -e $BINS_FILE ]; then
    echo "missing file: $BINS_FILE"
    echo "please run pipeline 'prepareBins' first"
    exit 1
fi

# count reads per sample in parallel by chromosome
zcat $BINS_FILE |
perl $MODULES_DIR/genome/sort_chroms.pl |
cut -f 1 |
uniq | # yields unique ordered chroms
parallel -j $N_CPU -k bash $ACTION_DIR/pivot/makePivot.sh {} |
awk 'BEGIN{OFS="\t"}{
    $4 = NR>1 ? (NR-1) : $4;
    print $0;
}' |
gzip -c > $COUNT_MATRIX_FILE
checkPipe

# report to log file
echo
echo $COUNT_MATRIX_FILE
echo "(columns 1-13 only)"
zcat $COUNT_MATRIX_FILE | cut -f 1-13 | head
echo "..."
zcat $COUNT_MATRIX_FILE | cut -f 1-13 | tail

# finish
echo "done"
echo
