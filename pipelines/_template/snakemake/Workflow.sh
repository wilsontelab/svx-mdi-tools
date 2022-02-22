#!/bin/bash

#--------------------------------------------------------------------
# sometimes your pipeline will require preparative steps
#--------------------------------------------------------------------

# demonstrate access to action scripts
source $ACTION_DIR/../do/echo.sh

# demonstrate use of a step module script
source $MODULES_DIR/demo/demo.sh

#--------------------------------------------------------------------
# snakemake actions typically defer all sequential steps to the Snakefile
#--------------------------------------------------------------------

# list files in a directory, make a plot, print a table, and sleep for 15 seconds

checkWorkflowStep 1 do-snakemake Snakefile

if [ "$SN_FORCEALL" != "" ]; then rm -rf $TASK_DIR/.snakemake; fi
snakemake $SN_DRY_RUN $SN_FORCEALL \
    --cores $N_CPU \
    --snakefile $ACTION_DIR/Snakefile \
    --directory $TASK_DIR \
    $DATA_NAME.directoryContents.txt $DATA_NAME.mtcars.csv
checkPipe

finishWorkflowStep

sleep 15
