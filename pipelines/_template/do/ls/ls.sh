#!/bin/bash

echo "listing directory contents..."

ls -l $INPUT_DIR > $DATA_FILE_PREFIX.directoryContents.txt
checkPipe
