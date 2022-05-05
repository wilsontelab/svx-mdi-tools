# action:
#     merge samples and add normalization data
# expects:
#     find2/find_svs.sh
# outputs:
#     $COVERAGE_PREFIX.rds

echo "assembling composite coverage file"
Rscript $ACTION_DIR/coverage/assemble_coverage.R
checkPipe

echo "done"
