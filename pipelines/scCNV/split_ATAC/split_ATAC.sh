
# extract the barcodes of good cells and number them
export GOOD_CELLS_FILE=$DATA_FILE_PREFIX.good_cells.csv
Rscript $ACTION_DIR/parse_good_cells.R

export BAM_HEADER=`samtools view -H $CELL_RANGER_BAM_FILE`

SORT_RAM_INT=$((TOTAL_RAM_INT - 4000000000))

samtools view $CELL_RANGER_BAM_FILE | 
perl $ACTION_DIR/filter_for_good_cells.pl | 
sort -k1,1n -k2,2 --compress-program=gzip --buffer-size=$SORT_RAM_INT"b" --parallel=$N_CPU | 
perl $ACTION_DIR/write_cell_bams.pl
checkPipe
