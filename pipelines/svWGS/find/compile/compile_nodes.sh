# action:
#     create node file indices used during SV event recontruction
# expects:
#     source $MODULES_DIR/scan/set_genome_vars.sh
#     source $MODULES_DIR/scan/set_alignment_vars.sh
#     extract/extract_nodes.sh
#     $MIN_CLIP
# input:
#     $EXTRACT_PREFIX.nodes.*.gz
# outputs:
#     $DATA_FILE_PREFIX.baseCoverage.*
#     insertSizes histogram

#-----------------------------------------------------------------
# define common actions
#-----------------------------------------------------------------
SORT_RAM_INT=`echo $TOTAL_RAM_INT | awk '{print int(($1 - 4000000000) / 2)}'`
SLURP_GZ="slurp -s 100M pigz -p $N_CPU -dc"
SORT="sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM_INT""b"
GROUP_BY="bedtools groupby"
PIGZ="pigz -p $N_CPU -c"
SLURP_OUT="slurp -s 100M -o"
SLURP_NODES="$SLURP_GZ $EXTRACT_PREFIX.nodes.*.gz"
PERL_COMPILE="perl $ACTION_DIR/compile"

#-----------------------------------------------------------------
# nodes columns
#-----------------------------------------------------------------
NODE=1 # node-level data
CLIP_LEN=2
CLIP_SEQ=3
NODE_CLASS=4
#---------------
JXN_TYPE=5 # edge/junction-level data
JXN_N=6
#---------------
FLAG=7 # alignment-level data
POS=8
MAPQ=9
CIGAR=10
SEQ=11
ALN_N=12
#---------------
MOL_ID=13 # molecule-level information  
MOL_CLASS=14
MOL_STRAND=15
IS_OUTER_CLIP1=16
IS_OUTER_CLIP2=17
#---------------
GAP=0 # node classes
SPLIT=1
OUTER_CLIP=2

#-----------------------------------------------------------------
echo "declaring junction edges"
#-----------------------------------------------------------------
function print_junctions {
    echo "  $2"
    $SLURP_NODES |
    awk '$'$JXN_N'>0&&$'$NODE_CLASS'=='$1'&&$'$JXN_TYPE'!="P"' |
    cut -f $NODE,$JXN_TYPE,$JXN_N,$MOL_ID |    
    $SORT -k4,4n -k3,3n -k1,1 |
    $GROUP_BY -g 3,4 -c 1,2 -o collapse,first |
    sed 's/,/\t/g' |

    $SORT -k1,1n -k2,2n -k3,3 -k4,4 -k5,5 | # JXN_N, MOL_ID, NODE1, NODE2, JXN_TYPE

    $PERL_COMPILE/set_junction_target_class.pl | # returns: node1, node2, type, targetClass, molId
    $SORT -k1,1 -k2,2 -k3,3 -k4,4 -k5,5n |


    $GROUP_BY -g 1,2,3,4 -c 5,5,5 -o collapse,count,count_distinct | # aggregate molId by node pair
    $PIGZ |
    $SLURP_OUT $COMPILE_PREFIX.$2.gz
    checkPipe
    echo "    done"
}
print_junctions $SPLIT 'sequenced_junctions'
print_junctions $GAP   'gap_junctions'

#-----------------------------------------------------------------
echo "indexing SV nodes by molecule"
#-----------------------------------------------------------------
$SLURP_NODES | # include all nodes, even unclipped outer
$SORT -k$MOL_ID,$MOL_ID"n" | 
$PERL_COMPILE/index_molecule_nodes.pl
checkPipe
echo "  done"

#-----------------------------------------------------------------
echo "indexing SV nodes by coordinate proximity"
#-----------------------------------------------------------------
$SLURP_NODES | # only include SV evidence nodes
awk 'BEGIN{OFS="\t"}$'$JXN_TYPE'!="*"||$'$CLIP_LEN'>='$MIN_CLIP'{
    split($'$NODE', x, ":");
    print x[1]":"x[2], x[3], $0;
}' |
$SORT -k1,1 -k2,2n | # i.e., CHROM_STRAND, POSITION
# $PIGZ |
# $SLURP_OUT $COMPILE_PREFIX.nodes_by_proximity.tmp.gz
# $SLURP_GZ $COMPILE_PREFIX.nodes_by_proximity.tmp.gz | 
$PERL_COMPILE/index_proximity.pl
checkPipe
echo "  done"

##-----------------------------------------------------------------
#echo "indexing SV nodes by node name"
##-----------------------------------------------------------------
#$SLURP_NODES |
#awk $AWK_SV_NODES | # only include SV evidence nodes
#$SORT -k$NODE,$NODE |
#$PERL_COMPILE/index_nodes.pl
#checkPipe

# clean up
rm -r $TMP_DIR_WRK

echo "done"
echo
