# action:
#     characterize novel SV paths in one or more input samples or nanopores
# input:
#     edges.rds file created by analyze_edges.R for one or more nanopore samples
# output:
#     data package

# characterize novel SV paths in one or more input samples or nanopores
mkdir -p $PLOTS_DIR
Rscript $FIND_MODULE_DIR/find_svs.R
checkPipe

echo "done"
