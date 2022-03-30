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

#-----------------------------------------------------------------
# define common actions
#-----------------------------------------------------------------
SORT_RAM_INT=`echo $TOTAL_RAM_INT | awk '{print int(($1 - 2000000000) / 2)}'`
SLURP_GZ="slurp -s 100M pigz -p $N_CPU -dc"
SORT="sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM_INT""b"
PERL_COMPILE="perl $ACTION_DIR/compile"
RSCRIPT_COMPILE="Rscript $ACTION_DIR/compile"

#-----------------------------------------------------------------
# node/junction columns
#-----------------------------------------------------------------
NODE_1=1 # node-level data
CLIP_LEN=2
CLIP_SEQ=3
#---------------
FLAG=4 # alignment-level data
POS=5
MAPQ=6
CIGAR=7
SEQ=8
ALN_N=9
#---------------
UMI=10
#---------------
NODE_CLASS=11
#---------------
JXN_TYPE=12 # edge/junction-level data
JXN_N=13
#---------------
MOL_ID=14 # molecule-level information  
IS_MERGED=15
IS_DUPLEX=16
STRAND_COUNT1=17
STRAND_COUNT2=18
MOL_CLASS=19
MOL_STRAND=20
IS_OUTER_CLIP1=21
IS_OUTER_CLIP2=22
TARGET_CLASS=23
SHARED_PROPER=24
#---------------
OUT_POS_1=25
OUT_POS_2=26
#---------------
SAMPLE=27
#---------------
NODE_2=28
# ...
#-----------------------------------------------------------------
JXN_NODES='$'$NODE_CLASS'!='$OUTER_CLIP'&&$'$JXN_TYPE'!="P"'

#-----------------------------------------------------------------
# determine if this is a single-sample find or a multi-sample compare
echo "setting find mode"
#-----------------------------------------------------------------
export FIND_MODE=find
export JUNCTION_FILES=`ls -1 $TASK_DIR/*compile.junction_edges.gz 2>/dev/null`
export LIBRARY_STAT_FILES=`ls -1 $TASK_DIR/*extract.library_stats.yml 2>/dev/null`
if [ "$JUNCTION_FILES" = "" ]; then
    export FIND_MODE=compare
    export JUNCTION_FILES=`ls -1 $TASK_DIR/*/*compile.junction_edges.gz 2>/dev/null`
    export LIBRARY_STAT_FILES=`ls -1 $TASK_DIR/*/*extract.library_stats.yml 2>/dev/null`
    if [ "$JUNCTION_FILES" = "" ]; then
        echo "could not find junction file(s) in:"
        echo "    $TASK_DIR"
        exit 1
    fi
fi
echo "  $FIND_MODE"

#-----------------------------------------------------------------
# parse the library names and MAX_TLEN values
echo "collecting library stats"
#-----------------------------------------------------------------
export SAMPLES=`cat $LIBRARY_STAT_FILES | grep SAMPLE | sed 's/SAMPLE:\s//'`
export MAX_TLENS=`cat $LIBRARY_STAT_FILES | grep MAX_TLEN | sed 's/MAX_TLEN:\s//'`

#-----------------------------------------------------------------
# begin to apply SV filters as possible
#-----------------------------------------------------------------
if [ "$ON_TARGET" = "2" ]; then 
    TARGET_CLASS_FILTER='$'$TARGET_CLASS'!~/-/'
elif [ "$ON_TARGET" = "1" ]; then
    TARGET_CLASS_FILTER='$'$TARGET_CLASS'!="--"'
else 
    TARGET_CLASS_FILTER=""
fi

#-----------------------------------------------------------------
# load files into shared memory for rapid access
#-----------------------------------------------------------------
source $MODULES_DIR/utilities/shell/create_shm_dir.sh
rm -f $SHM_DIR_WRK/*

# #   indexed nodes
# function load_nodes {
#     echo "loading $1 into RAM"
#     cp $COMPILE_PREFIX.$1.txt       $SHM_DIR_WRK/$1.txt
#     cp $COMPILE_PREFIX.$1.txt.index $SHM_DIR_WRK/$1.txt.index
# }
# # load_nodes nodes_by_proximity
# load_nodes outer_clips

#   genome
echo "loading $GENOME into RAM"
export SHM_GENOME_FASTA=$SHM_DIR_WRK/$GENOME.fa
cp $GENOME_FASTA     $SHM_GENOME_FASTA
cp $GENOME_FASTA.fai $SHM_GENOME_FASTA.fai

#-----------------------------------------------------------------
# break junctions in to continuity groups and parse to called SVs
echo "aggregating junction molecules into SV calls"
#-----------------------------------------------------------------
$SLURP_GZ $JUNCTION_FILES | 
awk 'BEGIN{OFS="\t"}'$TARGET_CLASS_FILTER'{
    split($'$NODE_1', n1, ":");
    split($'$NODE_2', n2, ":");
    print n2[1]":"n2[2]":"n1[1]":"n1[2], n1[3], n2[3], $0;
}' |
$SORT -k1,1 -k2,2n | 
perl $ACTION_DIR/find/group_junctions.pl | 
Rscript $ACTION_DIR/find/call_svs.R

exit 1

# #-----------------------------------------------------------------
# echo "indexing SV junction nodes by coordinate proximity" # used to collect junction-flanking evidence
# #-----------------------------------------------------------------
# $SLURP_NODES |
# awk 'BEGIN{OFS="\t"}'$JXN_NODES'{
#     split($'$NODE', x, ":");
#     print $0, x[1]":"x[2], x[3];
# }' |
# $SORT -k$PARTNER,$PARTNER -k$CHROM_STRAND,$CHROM_STRAND -k$POSITION,$POSITION"n" |  
# $MASK_NODES | 
# $PERL_COMPILE/index_proximity.pl
# checkPipe




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
