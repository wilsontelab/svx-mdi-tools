# action:
#     parse molecules, junctions, and paths from an assembled nodes file
# input:
#     processed nodes file, one line per alignment or junction edge
# output:
#     $ANALYZE_PREFIX.

# log file feedback
echo "analyzing molecule junctions and paths"

# # get ready
export NODES_FILE=$EXTRACT_PREFIX.nodes.txt.gz
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz
export SEQUENCES_FILE=$EXTRACT_PREFIX.sequences.txt.gz

# assemble the node path of every molecule, one line per segment
mkdir -p $PLOTS_DIR
Rscript $ACTION_DIR/analyze/analyze_nodes.R
checkPipe

echo "done"
