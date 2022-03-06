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
RSCRIPT_COMPILE="Rscript $ACTION_DIR/compile"
MASK_NODES="cat"

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
UMI=14
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
CHROM_STRAND=25
POSITION=26
#-----------------------------------------------------------------
# endpoint columns
#-----------------------------------------------------------------
EP_MOL_ID=1       # source molecule number
EP_MOL_CLASS=2    # P or V
EP_MOL_STRAND=3   # 0, 1 or 2 (for duplex)
EP_TARGET_CLASS=4 # TT tt t- etc.
EP_NODE=5 # umi:chromI:side:projPos:isSVClip
#-----------------------------------------------------------------
# node classes
#-----------------------------------------------------------------
GAP=0 # node classes
SPLIT=1
OUTER_CLIP=2

# some calculations only apply to collated molecules, i.e., svCapture
if [ "$IS_COLLATED" != "" ]; then
MASK_NODES="$PERL_COMPILE/mask_nodes.pl" # fill the SHARED_PROPER column using matchedProper.gz

#-----------------------------------------------------------------
echo "identifying SV molecules that share endpoints with proper molecules"
#-----------------------------------------------------------------
$SLURP_GZ $EXTRACT_PREFIX.endpoints.*.gz |
$SORT -k$EP_NODE,$EP_NODE | # sort by endpoint identity
$GROUP_BY -g $EP_NODE \
    -c $EP_MOL_ID,$EP_MOL_CLASS,$EP_MOL_CLASS \
    -o collapse,collapse,count_distinct | # aggregate molecules that claim each endpoint
awk '$NF>1&&$(NF-1)~/P/' | # more than one type of molecule claimed an endpoint id, and one was proper
perl -ne '
    chomp;
    my @f = split("\t");
    my @mId = split(",", $f[1]);
    my @mCl = split(",", $f[2]);
    map { $mCl[$_] ne "P" and print $mId[$_], "\n" } 0..$#mId
' | # retain all the variant molecules that matched a proper endpoint
$SORT -k1,1n |
$GROUP_BY -g 1 -c 1 -o count | # count the number of endpoints each SV molecule shared with a proper molecule
$PIGZ |
$SLURP_OUT $COMPILE_PREFIX.matchedProper.gz # save the result to act as a lookup for SV molecules to reject
checkPipe

#-----------------------------------------------------------------
echo "calculating duplex rates and strand family sizes"
#-----------------------------------------------------------------
cat \
<(cat $EXTRACT_PREFIX.strand_counts.1.txt | head -n1) \
<(cat $EXTRACT_PREFIX.strand_counts.*.txt | grep -v 'targetClass' | sort -k1,1 -k2,2 -k3,3n) |
$RSCRIPT_COMPILE/family_sizes.R
checkPipe

fi

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
$MASK_NODES |
$PERL_COMPILE/index_molecule_nodes.pl
checkPipe
echo "  done"

#-----------------------------------------------------------------
echo "indexing SV nodes by coordinate proximity"
#-----------------------------------------------------------------
$SLURP_NODES | # only include SV evidence nodes
awk 'BEGIN{OFS="\t"}$'$JXN_TYPE'!="*"||$'$CLIP_LEN'>='$MIN_CLIP'{
    split($'$NODE', x, ":");
    print $0, x[1]":"x[2], x[3];
}' |
$SORT -k$CHROM_STRAND,$CHROM_STRAND -k$POSITION,$POSITION"n" | 
# $PIGZ |
# $SLURP_OUT $COMPILE_PREFIX.nodes_by_proximity.tmp.gz
# $SLURP_GZ $COMPILE_PREFIX.nodes_by_proximity.tmp.gz | 
$MASK_NODES | 
$PERL_COMPILE/index_proximity.pl
checkPipe
echo "  done"


# clean up
rm -r $TMP_DIR_WRK

echo "done"
echo
