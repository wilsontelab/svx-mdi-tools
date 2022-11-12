# this script is run in parallel by CHROM
CHROM=$1

# extract good read on the query chromosome
slurp samtools view -q $MIN_MAPQ $COORDINATE_BAM_FILE $CHROM | # apply min quality filter
perl $ACTION_DIR/pivot/bam2bed.pl | # simple stream, will not change sort
bedtools intersect -v -a stdin -b <(zcat $BAD_REGIONS_FILE | awk '$1=="'$CHROM'"') | # apply bad region filter

# associate unique source molecules with large genome bins
bedtools intersect -sorted -wao -a <(zcat $BINS_FILE | awk '$1=="'$CHROM'"') -b stdin |

# create the cell pivot table
perl $ACTION_DIR/pivot/makePivot.pl $CHROM 
