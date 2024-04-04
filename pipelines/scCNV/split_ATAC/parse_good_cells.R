#=====================================================================================
# parse a table of good cell ids from singlecell.csv
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
suppressPackageStartupMessages({
    library(data.table)
})
#-------------------------------------------------------------------------------------
# parse environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
source(file.path(rUtilDir, 'utilities.R'))
checkEnvVars(list(
    string = c(
        'CELL_RANGER_CELL_FILE',
        'GOOD_CELLS_FILE'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2)
#=====================================================================================

#=====================================================================================
# extract good barcodes and number the cells
#-------------------------------------------------------------------------------------
x <- fread(env$CELL_RANGER_CELL_FILE)
isGoodCellCol <- "is__cell_barcode"
if(!(isGoodCellCol %in% names(x))) isGoodCellCol <- "is_cell_barcode"
isGoodCell <- x[[isGoodCellCol]] > 0
x <- x[isGoodCell][, .(barcode = barcode, cell_id = 1:.N)]
write.table(x, file = env$GOOD_CELLS_FILE, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = FALSE)
#=====================================================================================
