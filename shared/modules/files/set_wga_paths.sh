# action:
#     set file paths and prefixes for wgaSeq finding algorithms
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
# usage:
#     source $MODULES_DIR/files/set_wga_paths.sh

# wgaSeq do
export RATES_FILE=$DATA_GENOME_PREFIX.rates.txt
export COUNT_MATRIX_FILE=$DATA_GENOME_PREFIX.bins.size_$BIN_SIZE.bed.gz
export INS_SIZES_FILE=$DATA_GENOME_PREFIX.insertSizes.txt
export COUNTS_FILE=$DATA_GENOME_PREFIX.counts.txt
export SIMPLE_MANIFEST_FILE=$DATA_GENOME_PREFIX.manifest.csv
