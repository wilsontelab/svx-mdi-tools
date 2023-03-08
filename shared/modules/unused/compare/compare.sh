#!/bin/bash

#-----------------------------------------------------------------
# define common actions
#-----------------------------------------------------------------
SORT_RAM_INT=`echo $TOTAL_RAM_INT | awk '{printf "%0.0f", $1 - 2000000000}'`
SLURP_GZ="slurp -s 100M pigz -p $N_CPU -dc"
SORT="sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM_INT""b"
PERL_COMPARE="perl $ACTION_DIR"
RSCRIPT_COMPARE="Rscript $ACTION_DIR"

#-----------------------------------------------------------------
# SV columns
#-----------------------------------------------------------------
SV_ID=1
JUNCTION_NAME=2
MATCHING_NAMES=3
#-------------
TARGET_CLASS=4
JXN_TYPE=5
#-------------
CHROM_1=6
SIDE_1=7
POS_1=8
CHROM_2=9
SIDE_2=10
POS_2=11
#-------------
JXN_SEQ=12
MERGE_LEN=13
FAIDX_PADDING=14
GEN_REF_1=15
GEN_REF_2=16
#-------------
MICROHOM_LEN=17
MICROHOM_MATCH=18
JXN_BASES=19
SV_SIZE=20 
#-------------
N_TOTAL=21
N_SPLITS=22
N_GAPS=23
N_OUTER_CLIPS=24
N_DUPLEX=25
N_DUPLEX_ALL=26
NET_STRAND_COUNT=27
NET_STRAND_COUNT_ALL=28
N_SHARED_PROPER=29
N_SHARED_PROPER_ALL=30
#------------- 
IS_MERGED=31
SEQ_LEN=32
UMI=33
MAPQ=34
#-------------
CHUNK_OFFSET=35
CHUNK_SIZE=36
#-------------
SAMPLE=37

#-----------------------------------------------------------------
echo "getting largest insert size across all libraries"
#-----------------------------------------------------------------
export MAX_MAX_TLEN=`$PERL_COMPARE/getMaxMaxTLen.pl`

#-----------------------------------------------------------------
echo "comparing structural variants between samples"
#-----------------------------------------------------------------
$SLURP_GZ $FIND_GLOB_PREFIX.structural_variants.gz | 
$SORT -k$CHROM_2,$SIDE_2 -k$CHROM_1,$SIDE_1 -k$POS_1,$POS_1"n" | 
$PERL_COMPARE/makeGroups.pl | 
$RSCRIPT_COMPARE/compare.R
checkPipe

echo "done"
