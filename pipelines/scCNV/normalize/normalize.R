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
        'ACTION_DIR',
        'INPUT_DIR',
        'INPUT_NAME',
        'BAD_REGIONS_FILE',
        'PLOTS_DIR',
        'PLOT_PREFIX',
        'OUTPUT_FILE',
        'DATA_NAME',
        'MANIFEST_FILE' # created by us
    ),
    integer = c(
        "PLOIDY",
        "MAX_WINDOW_POWER",
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

# rep=  
# 3 = 24/LK-13, ~85%S, XX
# 70 = 69/LK-74, mid-S, XY + aneuploidy
# 37 = 32/LK-44, ?40%S, XY
# 8 = 79/LK-18, ~15%S, XX
# un = 
# 90 = 91/LK-92, 0%S, XX + gain chr22
# 81 = 81/LK-84, 0%S, XY + loss chr15

# cells <- mclapply(cell_ids, function(cell_id){
cells <- lapply("70", function(cell_id){ 
# cells <- lapply(as.character(c(3, 70, 37, 8, 90, 81)), function(cell_id){ 
    # message()
    # message(cell_id)
    cell <- fitCell(cell_id)
    # if(!is.null(cell$fit)) plotCellQC(cell_id, cell)
    # str(cell)
    
    cell
})

# dir <- "/nfs/turbo/umms-smithgd/pipelines/scripts/scCNV/human_embryo"
# file <- file.path(dir, "dev.rds")
# # saveRDS(P_rep_fs_gc, file = file)
# P_rep_fs_gc <- readRDS(file)

# plotReplicationProfiles()

stop("GOT TO HERE!")

# , mc.cores = env$N_CPU)
names(cells) <- cell_ids

# assemble and organize the per-cell data
message('recording cell metadata')
colData[, ':='(
    keep        =  cells[[cell_id]]$keep,
    rejected    = !cells[[cell_id]]$keep,
    windowPower = cells[[cell_id]]$windowPower,    
    sdLagDiff   = cells[[cell_id]]$sdLagDiff,
    replicating = cells[[cell_id]]$replicating,
    fractionS   = cells[[cell_id]]$fractionS,
    peakIsReplicated   = cells[[cell_id]]$replicationModel$peakIsReplicated,
    maxNR_unreplicated = cells[[cell_id]]$replicationModel$maxNR_unreplicated,
    modal_CN    = cells[[cell_id]]$modal_CN,
    ploidy      = cells[[cell_id]]$ploidy,
    ER_modal_CN = cells[[cell_id]]$ER_modal_CN,
    ER_ploidy   = cells[[cell_id]]$ER_ploidy,
    readsPerAllele = cells[[cell_id]]$readsPerAllele
), by = cell_id]
#=====================================================================================

#=====================================================================================
# save the interim output and manifest files
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
write.table(
    data.table(
        Project = env$DATA_NAME,
        Sample_ID = colData$cell_id,
        Description = if(is.null(colData$cell_name)) colData$cell_id else colData$cell_name,
        Keep = colData$keep
    ), 
    file = env$MANIFEST_FILE, 
    quote = FALSE, 
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)
#=====================================================================================
