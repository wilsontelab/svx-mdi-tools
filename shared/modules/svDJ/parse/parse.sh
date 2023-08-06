# action:
#     discover primer node positions
#     split chimeric amplicon fusions
#     purge duplex molecules
#     drop truncated reads that don't match primer ends
# input:
#     processed edge files, one line per alignment or junction edge
# output:
#     $PARSE_PREFIX.edges.rds

# score, filter and aggregate edges by molecule and junction
mkdir -p $PLOTS_DIR
Rscript $PARSE_MODULE_DIR/parse_amplicons.R
checkPipe

echo "done"
