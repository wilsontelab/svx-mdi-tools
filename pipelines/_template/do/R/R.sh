#!/bin/bash

echo "calling Rscript..."
Rscript $ACTION_DIR/R/my-script.R
checkPipe

echo "sleeping for 15 seconds"
sleep 15
