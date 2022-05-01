#!/bin/bash

# application-specific data_script called by workflow launcher
# organizes a set of serial pipeline steps, with staged execution

# validate genome file
checkWorkflowStep 1 checkBins _wgaSeq_align_main.sh

CHECK FOR FILES


 align reads to genome and create temporary output files
runWorkflowStep 2 align/align.sh

# make a table with cells counts per fixed-width bin
runWorkflowStep 3 count/count.sh

