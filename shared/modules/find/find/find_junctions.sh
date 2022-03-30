# action:
#     create node file indices used during SV event recontruction
# expects:
#     source $GENOMEX_MODULES_DIR/scan/set_genome_vars.sh
#     source $GENOMEX_MODULES_DIR/scan/set_alignment_vars.sh
#     extract/extract_nodes.sh
#     compile/compile_nodes.sh
# outputs:
#     $FIND_PREFIX.find.all_nodes.txt
#     $FIND_PREFIX.structural_variants.gz

echo "finding structural variants"

# load files into shared memory for rapid access
source $MODULES_DIR/utilities/shell/create_shm_dir.sh
rm -f $SHM_DIR_WRK/*

#   indexed nodes
function load_nodes {
    echo "loading $1 into RAM"
    cp $COMPILE_PREFIX.$1.txt       $SHM_DIR_WRK/$1.txt
    cp $COMPILE_PREFIX.$1.txt.index $SHM_DIR_WRK/$1.txt.index
}
load_nodes nodes_by_proximity
load_nodes outer_clips

#   genome
echo "loading $GENOME into RAM"
export SHM_GENOME_FASTA=$SHM_DIR_WRK/$GENOME.fa
cp $GENOME_FASTA     $SHM_GENOME_FASTA
cp $GENOME_FASTA.fai $SHM_GENOME_FASTA.fai

# filter molecules and find junctions
Rscript $ACTION_DIR/find/find_junctions.R
checkPipe

# add a mark to SV nodes that were also claimed by an earlier numbered SV
ALL_NODES_FILE=$FIND_PREFIX.all_nodes.txt
ALL_NODES_TMP=$ALL_NODES_FILE.tmp
mv -f $ALL_NODES_FILE $ALL_NODES_TMP
slurp -s 100M $ALL_NODES_TMP |
perl $ACTION_DIR/find/mark_molecule_repeats.pl |
slurp -s 100M -o $ALL_NODES_FILE
checkPipe

# clean up
rm -fr $SHM_DIR_WRK
rm -f  $ALL_NODES_TMP

echo "done"
