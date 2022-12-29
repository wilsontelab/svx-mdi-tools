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

# cell with shaping issues
# 39 = some hit, some miss on wavy "gains"
# 40 = NR_wm_reshaped - something is wrong with shape fit? loses peak? peak moves? FIXED!
# 59 = low coverage, passes, has a problem after shaping (prob.like 40) FIXED!
# 52 = extreme shape on chr16, exaggerates a gain? NO REAL ERROR, JUST THE VAGARIES OF SHAPE FITTING - suppress extreme shapes?

# cells with issues related to fractionS/replicating cell assignments
# 49 = very late replication (got)
# 51 = even later replication?? (missed)
# 93 = likely a miss of a late-replicating cell due to waviness, stays squashed
# 95 = possible very late replication, missed

# cells with replication segments misassignments (CN1_NAR1 instead of CN2_NAR0, etc.)
# others above
# 82 = late replicating with chr22+, some replication miscalls

# other miscellaneous cells
# 47 = very high res unreplicated
# 56 = small chr19 segmental, is it real?
# 83 = early replication with chr22+, nailed it
# 84 = similer to 83, somewhat later in S
# 86 = unreplicated (low res) partner to 83/84 (also 89, 90)
# 91 = low res, maybe early replicating, but it just gets squashed (note chr17, chr22), chr22+ still detectable

# the following cells are identified as replicating, with windowPower, cell$replicationModel$fractionS and final cell$fractionS from composite HMM
# 4       3       0.95    0.483793494006378 BAD?, probably a reject hyperseg, o/w need better early fit
# 6       3       0.7     0.666458326702586 GOOD
# 15      7       0.35    0.389833425515335 GOOD?, very low res, but probably right, re-analyze lower power value
# 3       2       0.825   0.755370942600468 GOOD
# 12      5       0.925   0.595459621837016 BAD?, reject hyperseg
# 29      5       0.125   0.110417206602506 GOOD
# 9       4       0.05    0.144906685705483 PROBABLY GOOD, or wrong due to waviness
# 8       3       0.25    0.260561737578124 GOOD
# 11      4       0.875   0.842922758227746 BAD?, reject hyperseg
# 17      4       0.9     0.781132451526303 PROBABLY GOOD, or wrong due to waviness
# 26      7       0.7     0.654049141142613 BAD, very low res, likely reject cell
# 37      3       0.55    0.544142647711136 GOOD
# 49      5       0.925   0.850760633879875 PROBABLY GOOD, a bit too low res, could re-analyze at power 4
# 70      4       0.7     0.668104049136526 GOOD
# 74      7       0.1     0.15496918560414  MAYBE GOOD, very low coverage cell
# 84      5       0.6     0.631453919762601 GOOD, best re-analyzed at lower power value (even 3)
# 92      7       0.05    0.039999596405234 PROBABLY GOOD, re-analyze at lower power value?
# 82      3       0.925   0.841817135150738 GOOD
# 83      4       0.5     0.526324580161477 GOOD

# TODO: plot NR frequency distribution for all window powers along with genome plot
#       many replicating cells are (not surprisingly) being pushed to too-high windowPower

# TODO: move to app, allow user overrides of key parameters, refit in app

# TODO: is there a suitable logic to implement post-HMM "correction" of replication-derived CNV segments?
# TODO: any way to automate badCell assignment when excessively wavy or hypersegmented?
# TODO: if the above are possible, then could go straight to common CNV/batch normalization in pipeline
#       if not, defer that to app after manual intervention

# TODO: establish how to set a replication timing value for bins/windows across many cells based on their fractionS + FAR data
#       see Excel in Gary Smith folder for the logic of the cell voting process
#       then use that to implement an external replication timing trainer as a replacement (adjunct?) to GC

# cell_ids <- as.character(c(3, 70, 37, 8, 90, 81))
# cell_ids <- cell_ids[1:2]
# cell_ids <- "92"
# cells <- lapply(cell_ids, fitCell)
cells <- mclapply(cell_ids, fitCell, mc.cores = env$N_CPU)
stop("GOT TO HERE!")

cells <- mclapply(cell_ids, fitCell, mc.cores = env$N_CPU)
names(cells) <- cell_ids

# assemble and organize the per-cell data
message('recording cell metadata')
getWindowsMetadata <- function(cell_id, key){
    windowPower <- cells[[cell_id]]$windowPower
    x <- cells[[cell_id]]$windows[[windowPower + 1]][[key]]
    if(is.null(x)) NA else x
}
getReplicationMetadata <- function(cell_id, key) {
    x <- cells[[cell_id]]$replicationModel[[key]]
    if(is.null(x)) NA else x
}
colData[, ':='(
    bad          =  cells[[cell_id]]$badCell,
    keep         = !cells[[cell_id]]$badCell,
    ploidy       = cells[[cell_id]]$ploidy,
    windowPower  = cells[[cell_id]]$windowPower,
    normLagDiffQ = getWindowsMetadata(cell_id, "normLagDiffQ"),     
    modal_NA     = cells[[cell_id]]$modal_NA,
    ER_modal_NA  = getWindowsMetadata(cell_id, "ER_modal_NA"),  
    ER_ploidy    = getWindowsMetadata(cell_id, "ER_ploidy"),  
    RPA          = getWindowsMetadata(cell_id, "RPA"), 
    replicating  = cells[[cell_id]]$cellIsReplicating,
    fractionS    = cells[[cell_id]]$fractionS,
    modelType    = getReplicationMetadata(cell_id, "modelType"),
    theta        = getReplicationMetadata(cell_id, "theta")
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
        WindowPower = colData$windowPower
    ), 
    file = env$MANIFEST_FILE, 
    quote = FALSE, 
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)
#=====================================================================================
