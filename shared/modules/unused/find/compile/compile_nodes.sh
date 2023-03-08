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
#     $DATA_FILE_PREFIX.junction_edges.*
#     $DATA_FILE_PREFIX.nodes_by_proximity.*
#     $DATA_FILE_PREFIX.outer_clips.*

#-----------------------------------------------------------------
# define common actions
#-----------------------------------------------------------------
SORT_RAM_INT=`echo $TOTAL_RAM_INT | awk '{printf "%0.0f", int(($1 - 4000000000) / 2)}'`
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
PARTNER=5 # chromI:side for the other node in a junction or molecule
#---------------
JXN_TYPE=6 # edge/junction-level data
JXN_N=7
#---------------
FLAG=8 # alignment-level data
POS=9
MAPQ=10
CIGAR=11
SEQ=12
ALN_N=13
#---------------
MOL_ID=14 # molecule-level information  
UMI=15
IS_MERGED=16
IS_DUPLEX=17
STRAND_COUNT1=18
STRAND_COUNT2=19
MOL_CLASS=20
MOL_STRAND=21
IS_OUTER_CLIP1=22
IS_OUTER_CLIP2=23
TARGET_CLASS=24
SHARED_PROPER=25
#---------------
OUT_POS_1=26
OUT_POS_2=27
#---------------
CHROM_STRAND=28
POSITION=29
#-----------------------------------------------------------------
# endpoint columns
#-----------------------------------------------------------------
EP_NODE=1      # umi:chromI:side:projPos:isSVClip
EP_MOL_ID=2    # source molecule number
EP_MOL_CLASS=3 # P or V
#-----------------------------------------------------------------
# node classes
#-----------------------------------------------------------------
GAP=0
SPLIT=1
OUTER_CLIP=2
#-----------------------------------------------------------------
# molecule classes
#-----------------------------------------------------------------
IS_PROPER='P'
IS_SV='V'
#-----------------------------------------------------------------
JXN_NODES='$'$NODE_CLASS'!='$OUTER_CLIP'&&$'$JXN_TYPE'!="P"'
CLIP_NODES='$'$NODE_CLASS'=='$OUTER_CLIP'&&$'$CLIP_LEN'>='$MIN_CLIP

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
<(zcat $EXTRACT_PREFIX.strand_counts.1.gz | head -n1) \
<(zcat $EXTRACT_PREFIX.strand_counts.*.gz | grep -v 'targetClass' | sort -k1,1 -k2,2 -k3,3n) |
$RSCRIPT_COMPILE/family_sizes.R
checkPipe

fi

#-----------------------------------------------------------------
echo "declaring SV junction edges"
#-----------------------------------------------------------------
$SLURP_NODES |
awk $JXN_NODES |
cut -f $NODE,$NODE_CLASS,$JXN_TYPE,$JXN_N,$MOL_ID |    
$SORT -k5,5n -k4,4n -k1,1 |
$GROUP_BY -g 4,5 -c 1,2,3 -o collapse,first,first | # group by molId+jxnN
sed 's/,/\t/g' |
$PERL_COMPILE/set_junction_target_class.pl | # returns: nodeClass, node1, node2, jxnType, targetClass, molId
$SORT -k1,1nr -k2,5 -k6,6n | # sort splits first
$GROUP_BY -g 1,2,3,4,5 -c 6,6,6 -o collapse,count,count_distinct | # aggregate molId by node pair
$PIGZ |
$SLURP_OUT $COMPILE_PREFIX.junction_edges.gz
checkPipe

#-----------------------------------------------------------------
echo "indexing SV junction nodes by coordinate proximity" # used to collect junction-flanking evidence
#-----------------------------------------------------------------
$SLURP_NODES |
awk 'BEGIN{OFS="\t"}'$JXN_NODES'{
    split($'$NODE', x, ":");
    print $0, x[1]":"x[2], x[3];
}' |
$SORT -k$PARTNER,$PARTNER -k$CHROM_STRAND,$CHROM_STRAND -k$POSITION,$POSITION"n" |  
$MASK_NODES | 
$PERL_COMPILE/index_proximity.pl
checkPipe

#-----------------------------------------------------------------
echo "indexing outer clip nodes by name" # used to add outer clips as SV evidence
#-----------------------------------------------------------------
$SLURP_NODES |
awk $CLIP_NODES |
$SORT -k$NODE,$NODE | 
$MASK_NODES | 
$PERL_COMPILE/index_clips.pl
checkPipe

# #-----------------------------------------------------------------
# echo "indexing outer clip nodes by molecule" # used to suppress molecule duplicates as unique SV evidence
# #-----------------------------------------------------------------
# $SLURP_NODES |
# awk $CLIP_NODES'&&$'$MOL_CLASS'=="'$IS_SV'"' |
# $SORT -k$MOL_ID,$MOL_ID"n" | 
# $PERL_COMPILE/index_clips.pl 'molecule'
# checkPipe

# #-----------------------------------------------------------------
# echo "indexing SV nodes by molecule"
# #-----------------------------------------------------------------
# $SLURP_NODES | # include all nodes in SV molecules, even unclipped outer
# $SORT -k$MOL_ID,$MOL_ID"n" | 
# $MASK_NODES |
# $PERL_COMPILE/index_molecule_nodes.pl
# checkPipe

# clean up
rm -r $TMP_DIR_WRK

echo "done"
echo
