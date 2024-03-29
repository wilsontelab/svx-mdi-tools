# action:
#     extract alignment information to create a coverage map and SV junction list
# expects:
#     source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/utilities/shell/create_temp_dir.sh
# input:
#     $ALIGNMENT_DIR
#     $ALIGNMENT_FILES
# outputs:
#     <bam/cram>.bins.gz # per cell
#     <bam/cram>.junctions.gz

# check for sufficent sort memory
if [ $RAM_PER_CPU_INT -lt 4000000000 ]; then
    echo
    echo "!!! the scCNV extract action requires --ram-per-cpu of at least 4G !!!"
    echo
    exit 1
fi

# run the alignments one sample at a time
echo "parsing read pairs from cell bam/cram files"
echo $ALIGNMENT_FILES | 
sed 's/ /\n/g' | 
parallel -j $N_CPU bash $ACTION_DIR/extract_cell.sh {}
checkPipe

echo "done"
