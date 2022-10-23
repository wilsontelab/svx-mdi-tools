#!/bin/bash

Rscript $ACTION_DIR/extract.R
checkPipe

echo "compressing QC plots to tar.gz"
cd $PLOTS_DIR
tar -czf $PLOTS_ARCHIVE *.qc.png
checkPipe
