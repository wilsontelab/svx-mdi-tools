#!/bin/bash

# Rscript $ACTION_DIR/TMP.R
Rscript $ACTION_DIR/normalize.R
checkPipe

# echo "compressing QC plots to tar.gz"
# cd $PLOTS_DIR/scCNV
# tar -czf $PLOTS_ARCHIVE *
# checkPipe
