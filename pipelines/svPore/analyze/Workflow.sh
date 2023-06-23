#!/bin/bash

# set derivative environment variables and file paths
export PAF_FILE_TYPE=target_reads
export GENOMEX_MODULES_DIR=$SUITES_DIR/genomex-mdi-tools/shared/modules
source $GENOMEX_MODULES_DIR/genome/set_genome_vars.sh
source $GENOMEX_MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/files/set_svx_paths.sh

# set the size units and thresholds
export USE_CHR_M=1

# set file paths, written by extract, read by analyze
export USAM_FILE=$DATA_FILE_PREFIX.$PAF_FILE_TYPE.unaligned.sam.gz
export EDGES_TMP_FILE=$REEXTRACT_PREFIX.edges.tmp.txt.gz # ~10K non-SV reads for for adapter training
export EDGES_SV_FILE=$REEXTRACT_PREFIX.edges.sv.txt.gz   # all SV reads

# pull SV reads and representative non-SV reads from pre-aligned long reads (high accuracy, CIGAR required)
export EXTRACT_PREFIX=$REEXTRACT_PREFIX
runWorkflowStep 1 extract $ACTION_DIR/../extract/reextract.sh

# perform read quality assessments and analyze novel SV paths
runWorkflowStep 2 analyze analyze.sh

rm -f $EDGES_TMP_FILE
