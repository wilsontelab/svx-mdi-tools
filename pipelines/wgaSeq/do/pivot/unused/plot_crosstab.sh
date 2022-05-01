#!/bin/bash

Rscript $PIPELINE_DIR/crosstab/plot_crosstab.R


#exit 666
#
#
## set files
#TMP_DIR_WRK=$TMP_DIR_LARGE/$PIPELINE_NAME"_crosstab_"$DATA_NAME
#mkdir -p $TMP_DIR_WRK
#rm -f $TMP_DIR_WRK/*
#TMP_FILE=$TMP_DIR_WRK/crosstab.tmp.bed
#TMP_HEADER=$TMP_FILE.header
#OUT_FILE="$OUTPUT_DIR/crosstab.$BIN_SIZE.csv"
#
## set common actions
#SLURP_GZ="slurp -s 100M pigz -p $N_CPU -dc"
#SORT="sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $MAX_SORT_RAM"
#PIGZ="pigz -p $N_CPU -c"
#SLURP_OUT="slurp -s 100M -o"
#INTERSECT="bedtools intersect -sorted"
#CROSSTAB_PERL="perl $PIPELINE_DIR/crosstab"
#AVERAGE="$CROSSTAB_PERL/groupBy.pl"
#
## create a composite exclusions file
#echo "creating composite exclusions file (gaps + bad regions)"
#EXLUSIONS_BED=$OUTPUT_DIR/crosstab.exlusions.bed.gz
#cat $GAP_FILE $BAD_REGIONS_FILE | bedutil collapse | $SORT -k1,1 -k2,2n | $PIGZ > $EXLUSIONS_BED
#checkPipe
#
## pull reads and filter against bad genome regions
## both possorted_bam.bam and BAD_REGIONS_FILE are sorted the same, as per samtools sort:
##   chr10..chr19,chr1,chr20..chr22,chr2..chr9,chrX,chrY
## this is not the same chrom order as linux sort -k1,1 -k2,2n
#echo "counting reads to create the bin x cell crosstab file"
#slurp -s 100M samtools view $CELL_RANGER_DIR/possorted_bam.bam |
#$CROSSTAB_PERL/bam2bed.pl | # simple stream, will not change sort
#$INTERSECT -c -a stdin -b $BAD_REGIONS_FILE |
#
## make a crosstab of chrom x bin as read pair counts
#$CROSSTAB_PERL/crosstab.pl > $TMP_FILE
#checkPipe
#
## determine the number of cells and set columns appropriately
#head -n1 $TMP_FILE > $TMP_HEADER
#N_CELLS=`cat $TMP_HEADER | awk '{print NF-6}'`
#MAX_DATA_COL0=$((6 + $N_CELLS))
#MAX_DATA_COL1=$(($MAX_DATA_COL0 + 1))
#MAX_DATA_COL2=$(($MAX_DATA_COL1 + 1))
#MAX_DATA_COL3=$(($MAX_DATA_COL2 + 1))
#BT_EXCL_COL=$(($MAX_DATA_COL0 + 6 + 1))
#BT_MAPP_COL=$(($MAX_DATA_COL1 + 5))
#BT_GC_COL=$(($MAX_DATA_COL2 + 5))
#
## determine the exlusions size, %GC and fraction mappability of each bin
## all bed streams are sorted lexicographically here to:
##   chr1,chr10..chr19,chr2,chr20..chr22,chr3..chr9,chrX,chrY
#echo "adding information about excluded bases, mappability and GC fraction"
#cat <(awk '{print $0"\texcluded_bases\tmappability\tgc_fraction"}' $TMP_HEADER) <(
#
#    # add exclusions info
#    tail -n+2 $TMP_FILE |
#    $SORT -k1,1 -k2,2n |
#    $INTERSECT -wao -a stdin -b $EXLUSIONS_BED |
#    cut -f 1-$MAX_DATA_COL0,$BT_EXCL_COL |
#    
#    # add mappability info within non-excluded bins
#    $INTERSECT -loj -a stdin -b <(
#        $SLURP_GZ $MAPPABILITY_FILE |
#        $SORT -k1,1 -k2,2n |
#        $INTERSECT -v -a stdin -b $EXLUSIONS_BED
#    ) |
#    $AVERAGE $MAX_DATA_COL1 $BT_MAPP_COL |
#    
#    # add GC percent info within non-excluded bins
#    $INTERSECT -loj -a stdin -b <(
#        $SLURP_GZ $GC_FILE |
#        $SORT -k1,1 -k2,2n |
#        $INTERSECT -v -a stdin -b $EXLUSIONS_BED
#    ) |
#    $AVERAGE $MAX_DATA_COL2 $BT_GC_COL
#) |
#sed 's/\t/,/g' > $OUT_FILE
#checkPipe
#
## clean up
#rm -r $TMP_DIR_WRK

echo "done"
echo


#/treehouse/wilsonte_lab/umms-glover/data/linCNV/projects/GarySmith/GarySmith_100720/exlusions.bed.gz
#/treehouse/wilsonte_lab/ssd/genomes/Blacklist/lists/hg38-blacklist.v2.bed
#/home/wilsonte_lab/clubhouse/genomes/hg38/hg38.gap.bed
#/home/wilsonte_lab/clubhouse/genomes/hg38/hg38.kmer_50.bin_1000.bed.gz
#/home/wilsonte_lab/clubhouse/genomes/hg38/hg38.gc5Base.bin_1000.bed.gz

