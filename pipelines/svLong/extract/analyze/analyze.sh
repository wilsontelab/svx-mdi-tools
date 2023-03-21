# action:
#     parse molecules, junctions, and paths from an assembled nodes file
# input:
#     processed nodes file, one line per alignment or junction edge
# output:
#     $ANALYZE_PREFIX.

# log file feedback
echo "analyzing molecules junctions and paths"

# # get ready
export NODES_FILE=$EXTRACT_PREFIX.nodes.txt.gz
export COVERAGE_FILE=$EXTRACT_PREFIX.windowCoverage.txt.bgz
# SORT_RAM=$(( $MAX_SORT_RAM_INT / 2 )) 
# SORT_RAM=$SORT_RAM"b"
# SORT="sort --parallel $N_CPU -T $TMP_DIR_WRK -S $SORT_RAM --compress-program=pigz"

# assemble the node path of every molecule, one line per segment
Rscript $ACTION_DIR/analyze/analyze_nodes.R
