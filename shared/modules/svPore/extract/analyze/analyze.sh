# action:
#     parse molecules, junctions, and paths from assembled edges files
# input:
#     processed edges edges files, one line per alignment or junction edge
# output:
#     $ANALYZE_PREFIX.edges.rds

# score, filter and aggregate edges by molecule and junction
mkdir -p $PLOTS_DIR
Rscript $EXTRACT_STEP_DIR/analyze_edges.R
checkPipe

echo "done"
