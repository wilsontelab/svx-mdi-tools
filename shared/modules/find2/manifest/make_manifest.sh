# action:
#     create a sample manifest for use by the svx family of apps
# expects:
#     find2/find_svs.sh
# outputs:
#     $MANIFEST_PREFIX.distribution_plots.zip
#     $MANIFEST_PREFIX.sample_manifest.csv

echo "zipping distribution plots"
ZIP_FILE=$MANIFEST_PREFIX.distribution_plots.zip
if [ "$FIND_MODE" = "find" ]; then
    PLOT_GLOB=$TASK_DIR/plots/*.png
else
    PLOT_GLOB=$TASK_DIR/*/plots/*.png
fi
zip -j $ZIP_FILE $PLOT_GLOB
checkPipe

echo "assembling sample manifest table"
Rscript $ACTION_DIR/manifest/make_manifest.R
checkPipe

echo "done"
