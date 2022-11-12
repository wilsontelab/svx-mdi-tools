# convert one cell's bam alignments to coverage map and SV junction files

# set operating parameters
ALIGNMENT_FILE=$1
SORT_RAM=$(( $RAM_PER_CPU_INT - 3000000000 ))
SORT_RAM=$(( $SORT_RAM / 4 ))
SORT_RAM=$SORT_RAM"b"

# extract bin and junction data
perl $ACTION_DIR/extract.pl $ALIGNMENT_FILE | 

# purge duplicate read pairs
sort -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz -k1,1 -k2,2n | 
bedtools groupby -g 1 -c 3,4 -o first,first | 

# make a file of SV junction edges
# each junction present twice for indexing purposes
tee >(
    cut -f 3 | 
    sed -e 's/:::/\n/g' -e 's/::/\t/g' -e 's/:/\t/g' |
    awk 'BEGIN{OFS="\t"; x=1}$1!="X"&&$1!=0&&$4!=0{ print $0, x; print $4, $5, $6, $1, $2, $3, $7, x; x++}' |
    sort -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz -k1,1n -k3,3n -k4,4n -k6,6n |
    gzip -c > $ALIGNMENT_FILE.junctions.gz
) | 

# make a file of genomic bin counts
cut -f2 | 
sed -e 's/::/\n/g' -e 's/:/\t/g' |
awk '$1!=0' | 
sort -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz -k1,1n -k2,2n |
bedtools groupby -g 1,2 -c 2 -o count |
gzip -c > $ALIGNMENT_FILE.bins.gz
