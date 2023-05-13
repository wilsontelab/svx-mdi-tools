# action:
#     parse molecules, junctions, and paths from an assembled edges file
# input:
#     processed edges file, one line per alignment or junction edge
# output:
#     $ANALYZE_PREFIX.edges.rds

# score, filter and aggregate edges by molecule and junction
mkdir -p $PLOTS_DIR
Rscript $ACTION_DIR/analyze/analyze_edges.R
checkPipe

echo "done"
