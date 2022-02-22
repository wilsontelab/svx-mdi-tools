#!/bin/bash

#--------------------------------------------------------------------
# sometimes your pipeline will require preparative steps
#--------------------------------------------------------------------

# demonstrate access to action scripts
source $ACTION_DIR/echo.sh

# demonstrate use of a step module script
source $MODULES_DIR/demo/demo.sh

#--------------------------------------------------------------------
# most of a pipeline is typically coordinated into ordered steps
#--------------------------------------------------------------------

# list files in a directory
runWorkflowStep 1 do-ls ls/ls.sh

# make a plot, print a table, and sleep for 15 seconds
runWorkflowStep 2 do-R R/R.sh
