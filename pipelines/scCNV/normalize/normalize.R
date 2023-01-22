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
        'MANIFEST_FILE_IN', # optional user-provided manifest
        'MANIFEST_FILE',    # created by us
        'SHAPE_CORRECTION'        
    ),
    integer = c(
        "PLOIDY",
        "MIN_WINDOW_POWER",
        "MAX_WINDOW_POWER",
        "MAX_FIT_WINDOWS",
        "N_CPU"
    ),
    double = c(
        "N_SD_HALFCN",
        "MAX_NA95",
        "MIN_MAPPABILITY",
        "MIN_FRACTION_S",
        "KEEP_THRESHOLD"
    ),
    logical = c(
        "VERBOSE_PLOTS"
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
sourceScripts(scCnvSharedDir, c('loadHdf5', 'parseGenomeBins', 'fitCells', 'normalizeBatch'))
source(file.path(rUtilDir, 'smooth.R'))
#=====================================================================================

#=====================================================================================
# characterize the cells
#-------------------------------------------------------------------------------------
message('characterizing individual cells')
# cell_ids <- as.character(c(53,126,145,317,298,316,79,141,144,16))
# cell_ids <- as.character(c(3, 70, 37, 8, 90, 81, 78))
# cell_ids <- "0"
# cells <- lapply(cell_ids, fitCell_1)
cells <- mclapply(cell_ids, fitCell_1, mc.cores = env$N_CPU)
names(cells) <- cell_ids

# stop("XXXXXXXXXXXXXXXXX")
DIR <- env$TASK_ACTION_DIR
RDS_FILE <- file.path(DIR, "TEST.rds")
saveRDS(cells, file = RDS_FILE)
# cells <- readRDS(RDS_FILE)
# stop("XXXXXXXXXXXXXXXXX")

# assemble and organize the per-cell data
message('recording cell metadata')
shapeKey <- if(env$SHAPE_CORRECTION %in% c('cell', 'both')) "shaped" else "unshaped"
getCellValue <- function(cell_id, key){
    x <- cells[[cell_id]][[key]]
    if(is.null(x)) NA else x 
}
getWindowsMetadata <- function(cell_id, key){
    if(cells[[cell_id]]$badCell) return(NA)
    x <- cells[[cell_id]]$windows[[shapeKey]][[key]]
    if(is.null(x)) NA else x
}
getReplicationMetadata <- function(cell_id, key) {
    if(cells[[cell_id]]$badCell) return(NA)
    x <- cells[[cell_id]]$replicationModel[[shapeKey]][[key]]
    if(is.null(x)) NA else x
}
colData[, ':='( # initialize data types (in case first cell is NA)
    #---------------------------------- cell status information
    bad          = TRUE, # cells rejected BEFORE detailed analysis, no fit or model is present
    keep         = TRUE, # cells marked as questionable AFTER detailed analysis
    #---------------------------------- coverage and variance details
    ploidy       = 0L,
    windowPower  = 0L,
    modal_NA     = 0L,
    nrModelType  = "NA",  
    NA95         = 0.0,  
    ER_modal_NA  = 0.0,  
    ER_ploidy    = 0.0,  
    RPA          = 0.0,
    cnsd         = 0.0,
    #---------------------------------- replication model details
    replicating  = TRUE,
    fractionS    = 0.0,
    repGcRatio   = 0.0,
    modelType    = "NA",
    theta        = 0.0
)]
updateColData <- function(){
    colData[, ':='(
        #---------------------------------- cell status information
        bad          = getCellValue(cell_id, "badCell"), # cells rejected BEFORE detailed analysis, not fit or model is present
        keep         = getCellValue(cell_id, "keep"),    # cells marked as questionable AFTER detailed analysis
        #---------------------------------- coverage and variance details
        ploidy       = getCellValue(cell_id, "ploidy"),
        windowPower  = getCellValue(cell_id, "windowPower"),
        modal_NA     = getCellValue(cell_id, "modal_NA"),
        nrModelType  = getWindowsMetadata(cell_id, "nrModelType"),    
        NA95         = getWindowsMetadata(cell_id, "NA95"),  
        ER_modal_NA  = getWindowsMetadata(cell_id, "ER_modal_NA"),  
        ER_ploidy    = getWindowsMetadata(cell_id, "ER_ploidy"),  
        RPA          = getWindowsMetadata(cell_id, "RPA"),
        cnsd         = getCellValue(cell_id, "cnsd"),
        #---------------------------------- replication model details
        replicating  = getCellValue(cell_id, "cellIsReplicating"),
        fractionS    = getCellValue(cell_id, "fractionS"),
        repGcRatio   = getReplicationMetadata(cell_id, "gcRatio"),
        modelType    = getReplicationMetadata(cell_id, "modelType"),
        theta        = getReplicationMetadata(cell_id, "theta")
    ), by = cell_id]    
}
updateColData()
#=====================================================================================

#=====================================================================================
# parse the sample manifest; one run of normalize.R applies to a single Project, 
# which might have multiple Samples, each with multiple Cells
# if manifest not provided, all Cells are assumed to derive from the same Sample
#-------------------------------------------------------------------------------------
manifest <- if(env$MANIFEST_FILE_IN != "NA"){
    if(!file.exists(env$MANIFEST_FILE_IN)) stop(paste("unknown manifest file:", env$MANIFEST_FILE_IN))
    m <- fread(env$MANIFEST_FILE_IN)
    if(!all(c("Sample_Name","Cell_Name") %in% names(m))) stop("bad manifest file; expected columns Sample_Name and Cell_Name")
    if(anyDuplicated(m$Cell_Name)) stop("bad manifest file; all Cell_Name entries must be unique")
    if(is.null(colData$cell_name)) stop("don't know how to interpret manifest$Cell_Name because colData$cell_name is NULL")
    if(!all(m$Cell_Name %in% colData$cell_name)) stop("bad manifest file; one or more Cell_Name missing from colData$cell_name")
    if(!all(colData$cell_name %in% m$Cell_Name)) stop("bad manifest file; one or more colData$cell_name values are missing from manifest")
    setkey(m, "Cell_Name")
    names(m)[names(m) %in% c("", "", "", "")]
    x <- data.table(
        Project     = env$DATA_NAME,
        Sample_ID   = colData$cell_id, # the index value, i.e. "0", "1", etc.
        Sample_Name = m[colData$cell_name, Sample_Name],
        Cell_Name   = m[colData$cell_name, Cell_Name],
        Description = m[colData$cell_name, if(is.null(m$Description)) paste(Sample_Name, Cell_Name, sep = ":") else Description]
    )
    cbind(x, m[colData$cell_name, .SD, .SDcols = names(m)[!(names(m) %in% names(x))]])
} else {
    cellNames <- if(is.null(colData$cell_name)) colData$cell_id else colData$cell_name
    data.table(
        Project     = env$DATA_NAME,
        Sample_ID   = colData$cell_id,
        Sample_Name = env$DATA_NAME,
        Cell_Name   = cellNames,
        Description = cellNames
    ) 
}
manifest[, ":="(
    Yield = colData$total_num_reads,
    Quality = colData$windowPower + colData$bad
)]
#=====================================================================================

#=====================================================================================
# compare cells across a batch, i.e., a sample
#-------------------------------------------------------------------------------------
sampleNames <- unique(manifest$Sample_Name)
if(env$SHAPE_CORRECTION %in% c('sample', 'both')){ # normalize windows for batch effects
    message('calculating multi-cell batch correction, i.e., common shape by sample')
    for(sampleName in sampleNames){
        cell_ids <- manifest[Sample_Name == sampleName, Sample_ID]
        cells[cell_ids] <- normalizeBatch(sampleName, cells[cell_ids])
    }

    message('applying multi-cell batch correction to individual cells')
    cells <- mclapply(cells, fitCell_2, mc.cores = env$N_CPU) 
    shapeKey <- "batched"
    updateColData()
}
sampleProfiles <- lapply(sampleNames, function(sampleNames){ # determine the sample's reference CN
    cell_ids <- manifest[Sample_Name == sampleNames, Sample_ID]
    getSampleProfile(cells[cell_ids])
})
names(sampleProfiles) <- sampleNames
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
    raw_counts = raw_counts,
    states = list(
        replicating = replicatingStates,
        nonReplicating = nonReplicatingStates
    ),
    sampleProfiles = sampleProfiles
), file  = env$OUTPUT_FILE)
write.table(
    manifest, 
    file = env$MANIFEST_FILE, 
    quote = TRUE, 
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)
#=====================================================================================
