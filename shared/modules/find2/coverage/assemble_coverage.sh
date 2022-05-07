# action:
#     merge samples and add normalization data
# expects:
#     find2/find_svs.sh
# outputs:
#     $COVERAGE_PREFIX.rds

# coverage map only applies to whole genomic SV analysis
if [[ "$TARGETS_BED" = "" || "$TARGETS_BED" = "NA" || "$TARGETS_BED" = "null" ]]; then

echo "assembling composite coverage file"
Rscript $ACTION_DIR/coverage/assemble_coverage.R
checkPipe

echo "done"

fi
