# action:
#     create node file indices used during SV event recontruction
# expects:
#     source $MODULES_DIR/scan/set_genome_vars.sh
#     source $MODULES_DIR/scan/set_alignment_vars.sh
#     extract/extract_nodes.sh
#     compile/compile_nodes.sh
# input:
#     
# outputs:
#     

echo "finding structural variants"

# load files into shared memory for rapid access
# source $MODULES_DIR/utilities/shell/create_shm_dir.sh

#########
source $MODULES_DIR/utilities/shell/create_temp_dir.sh
export SHM_DIR_WRK=$TMP_DIR_WRK

rm -f $SHM_DIR_WRK/*
#   indexed nodes
function load_nodes {
    echo "loading $1 into RAM"
    cp $COMPILE_PREFIX.$1.txt       $SHM_DIR_WRK/$1.txt
    cp $COMPILE_PREFIX.$1.txt.index $SHM_DIR_WRK/$1.txt.index
}
load_nodes nodes_by_molecule
load_nodes nodes_by_proximity
#   genome
echo "loading $GENOME into RAM"
export SHM_GENOME_FASTA=$SHM_DIR_WRK/$GENOME.fa
cp $GENOME_FASTA     $SHM_GENOME_FASTA
cp $GENOME_FASTA.fai $SHM_GENOME_FASTA.fai

########################
echo $SHM_DIR_WRK
ls -lh $SHM_DIR_WRK
exit 1



# filter molecules and find junctions
Rscript $ACTION_DIR/find/find_junctions.R
checkPipe

# add a mark to SV nodes that were also claimed by an earlier numbered SV
ALL_NODES_FILE=$FIND_PREFIX.all_nodes.txt
mv -f $ALL_NODES_FILE $ALL_NODES_FILE.tmp
slurp -s 100M $ALL_NODES_FILE.tmp |
perl $ACTION_DIR/find/mark_molecule_repeats.pl |
slurp -s 100M -o $ALL_NODES_FILE
checkPipe

# clean up
rm -fr $SHM_DIR_WRK
rm -f $ALL_NODES_FILE.tmp

echo "done"
