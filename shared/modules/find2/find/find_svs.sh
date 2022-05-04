# action:
#     create node file indices used during SV event recontruction
# expects:
#     source $GENOMEX_MODULES_DIR/scan/set_genome_vars.sh
#     extract/extract_nodes.sh
#     compile/compile_nodes.sh
# outputs:
#     $FIND_PREFIX.metadata.yml
#     $FIND_PREFIX.target_regions.bed
#     $FIND_PREFIX.structural_variants.gz
#     $FIND_PREFIX.junction_molecules.gz

#-----------------------------------------------------------------
# define common actions
#-----------------------------------------------------------------
SLURP_GZ="slurp -s 100M pigz -p $N_CPU -dc"
SORT="sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM_INT""b"

#-----------------------------------------------------------------
# node/junction columns
#-----------------------------------------------------------------
NODE_1=1 # node-level data
CLIP_LEN_1=2
CLIP_SEQ_1=3
#---------------
FLAG_1=4 # alignment-level data
POS_1=5
MAPQ_1=6
CIGAR_1=7
SEQ_1=8
ALN_N_1=9
#---------------
UMI_1=10
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
#===============
NODE_2=28
# ...

#-----------------------------------------------------------------
# parse the library names and MAX_TLEN values
echo "collecting library stats"
#-----------------------------------------------------------------
export SAMPLES=`cat $LIBRARY_STAT_FILES | grep SAMPLE | sed 's/SAMPLE:\s//'`
export MAX_TLENS=`cat $LIBRARY_STAT_FILES | grep MAX_TLEN | sed 's/MAX_TLEN:\s//'`
export SAMPLES=`echo "$SAMPLES" | tr '\n' ' '`
export MAX_TLENS=`echo "$MAX_TLENS" | tr '\n' ' '`
echo "SAMPLES: $SAMPLES" > $FIND_PREFIX.metadata.yml
echo "MAX_TLENS: $MAX_TLENS" >> $FIND_PREFIX.metadata.yml

#-----------------------------------------------------------------
# begin to apply SV filters
#-----------------------------------------------------------------
if [ "$ON_TARGET" = "2" ]; then 
    TARGET_CLASS_FILTER='$'$TARGET_CLASS'!~/-/'
elif [ "$ON_TARGET" = "1" ]; then
    TARGET_CLASS_FILTER='$'$TARGET_CLASS'!="--"'
else 
    TARGET_CLASS_FILTER=""
fi

#-----------------------------------------------------------------
# break junctions into continuity groups and parse to called SVs
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
checkPipe

# clean up
rm -fr $SHM_DIR_WRK

echo "done"
