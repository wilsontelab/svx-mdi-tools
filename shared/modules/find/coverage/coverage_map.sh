# action:
#     create a coverage map of the genome and plot DNA fragment lengths
# expects:
#     source $MODULES_DIR/scan/set_genome_vars.sh
#     source $MODULES_DIR/scan/set_alignment_vars.sh
#     extract/extract_nodes.sh
# input:
#     $EXTRACT_PREFIX.spans.*.gz
# outputs:
#     $DATA_FILE_PREFIX.baseCoverage.*
#     insertSizes histogram

# coverage map only applies to whole genomic SV analysis
if [ "$TARGETS_BED" = "" ]; then

# set the sort parameters
SORT_RAM_INT=`echo $TOTAL_RAM_INT | awk '{print int(($1 - 4000000000) / 2)}'`

# extract SV info from name sorted bam and send all reads to coordinate sort
echo "processing genome fragment spans into coverage map"
slurp -s 500M pigz -p $N_CPU -dc $EXTRACT_PREFIX.spans.*.gz |
sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM_INT"b" | # the de-duplication sort
bedtools groupby -g 1,2,3,4 -c 5 -o first | # the de-duplication grouping
perl $ACTION_DIR/coverage/parse_spans.pl | # convert to 3-column BED + type in col 4
sort --parallel=$N_CPU -T $TMP_DIR_WRK -S $SORT_RAM_INT"b" -k1,1n -k2,2n | # the span coordinate sort
# pigz -p $N_CPU -c | 
# slurp -s 500M -o $DATA_FILE_PREFIX.baseCoverage.tmp.gz
# slurp -s 500M pigz -p $N_CPU -dc $DATA_FILE_PREFIX.baseCoverage.tmp.gz | 
perl $ACTION_DIR/coverage/base_coverage.pl
checkPipe

# clean up
rm -rf $TMP_DIR_WRK/*
# rm -rf $EXTRACT_PREFIX.spans.*.gz

echo "done"
echo

fi
