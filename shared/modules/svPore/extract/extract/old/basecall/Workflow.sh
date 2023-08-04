#!/bin/bash

# set derivative environment variables
export SHARED_SUITE=genomex-mdi-tools
export SHARED_MODULES_DIR=$SUITES_DIR/$SHARED_SUITE/shared/modules
export SHARED_MODULE_DIR=$SHARED_MODULES_DIR/nanopore/basecall

# set the dorado directories
echo "check nanopore directories"
DORADO_DIR=${NANOPORE_DIR}/dorado
DORADO_VERSION_NAME=dorado-${DORADO_VERSION}
DORADO_VERSION_DIR=${DORADO_DIR}/${DORADO_VERSION_NAME}
DORADO_EXECUTABLE=${DORADO_VERSION_DIR}/bin/dorado
if [ ! -f ${DORADO_EXECUTABLE} ]; then
    echo "missing Dorado executable: "${DORADO_EXECUTABLE}
    exit 1
fi 

# set ONT model paths
ONT_MODELS_DIR=${NANOPORE_DIR}/models
ONT_MODEL_DIR=${ONT_MODELS_DIR}/${ONT_MODEL}
if [ ! -d ${ONT_MODEL_DIR} ]; then
    echo "missing ONT model: "${ONT_MODEL_DIR}
    exit 1
fi 

# convert ONT read files from POD5 to FASTQ, i.e., call bases
runWorkflowStep 1 basecall $SHARED_MODULE_DIR/basecall.sh
