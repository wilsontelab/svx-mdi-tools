# extract bin count data from the Cell Ranger DNA hdf5 output

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(data.table)
library(rhdf5)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'INPUT_DIR',
        'INPUT_NAME',
        'TASK_DIR',
        'OUTPUT_DIR',
        'DATA_NAME'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#=====================================================================================


#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
h5File <- file.path(env$INPUT_DIR, env$INPUT_NAME, "cnv_data.h5")
if(!file.exists(h5File)) stop(paste("file not found:", h5File)
h5 = H5Fopen(h5File)
str(h5)
print(h5)
#=====================================================================================

stop("XXXXXXXXXXXXXXXXXXXX")
