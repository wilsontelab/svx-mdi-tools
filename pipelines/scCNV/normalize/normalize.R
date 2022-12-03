#=====================================================================================
# reformat data from step 'extract', or 10x Cell Ranger DNA, as input
# analyze cells one at a time to set windows, fit a GC bias model, and run a copy number HMM
#-------------------------------------------------------------------------------------
# regardless of the input data source, data are in HDF5 format as per 10x Cell Ranger DNA
#    https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/output/hdf5
#    also see accompanying summary file cell_ranger_file_format.txt
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
message("initializing")
suppressPackageStartupMessages({
    library(rhdf5)
    library(data.table)
    library(parallel)
})
#-------------------------------------------------------------------------------------
# parse environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
source(file.path(rUtilDir, 'utilities.R'))
checkEnvVars(list(
    string = c(
        'INPUT_DIR',
        'INPUT_NAME',
        'BAD_REGIONS_FILE',
        'PLOTS_DIR',
        'PLOT_PREFIX',
        'OUTPUT_FILE'        
    ),
    integer = c(
        "PLOIDY",
        "MAX_WINDOW_BINS",
        "N_CPU"
    ),
    double = c(
        "N_SD_HALFCN",
        "MIN_MAPPABILITY",
        "TRANSITION_PROBABILITY"
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
sourceClass(env, 'nbinomCountsGC2')
sourceClass(env, 'hmmEPTable')
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2)
#-------------------------------------------------------------------------------------
# set some constants and paths
dir.create(env$PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
unlink(paste(env$PLOT_PREFIX, "*", sep = "."))
#=====================================================================================

#=====================================================================================
# perform standard actions to load and parse genome bins
#-------------------------------------------------------------------------------------
scCnvSharedDir <- file.path(env$MODULES_DIR, 'scCNV')
sourceScripts(scCnvSharedDir, c('loadHdf5', 'parseGenomeBins', 'fitCells'))
#=====================================================================================

#=====================================================================================
# characterize the individual cells
# first pass fit at the bin resolution expected for Poisson without over/under-dispersion
#-------------------------------------------------------------------------------------
message('characterizing individual cells')
colData[, modal_CN := env$PLOIDY]
cells <- mclapply(cell_ids, function(cell_id){
    # message()
    # message(cell_id)
    cell <- fitCell(cell_id)
    if(!is.null(cell$fit)) plotCellQC(cell_id, cell)
    cell
}, mc.cores = env$N_CPU)
names(cells) <- cell_ids

# assemble and organize the per-cell data
message('recording cell metadata')
colData[, window_size := 0L] # needed to prevent errors if first cell has window_size == NA
colData[, ':='(
    rejected    = cells[[cell_id]]$rejected,
    stage       = cells[[cell_id]]$stage,
    pass        = cells[[cell_id]]$pass,
    window_size = cells[[cell_id]]$window_size
), by = cell_id]
#=====================================================================================

#=====================================================================================
# save the interim output file
#-------------------------------------------------------------------------------------
saveRDS(list(
    env = env,
    constants = constants,
    metadata = metadata,
    rowRanges = rowRanges, 
    windows = windows,
    colData = colData,
    cells = cells,
    raw_counts = raw_counts
), file  = env$OUTPUT_FILE)
#=====================================================================================
