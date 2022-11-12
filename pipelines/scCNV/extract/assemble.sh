# action:
#     combine all cells into single composite files for the entire sample
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $ALIGNMENT_DIR
#     $ALIGNMENT_FILES
# outputs:
#     $DATA_GENOME_PREFIX.extract.h5 # same format as 10x scCNV pipeline cnv_data.h5
#     $DATA_GENOME_PREFIX.extract.junctions.gz

# assemble bin coverage
echo "making 10x scCNV-compatible hdf5 file with bin coverage data"
Rscript $ACTION_DIR/assemble_bins.R

# TODO: assemble junctions

echo "done"
