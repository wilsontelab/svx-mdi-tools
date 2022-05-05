
# merge samples and add normalization data

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(data.table)
library(yaml)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOMES_DIR',
        'GENOME',
        'FIND_PREFIX',
        'FIND_MODE',
        'TASK_DIR',
        'OUTPUT_DIR',
        'DATA_NAME'
    ), 
    integer = c(
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# parse the project name and directory
if(env$FIND_MODE == "find"){
    projectName <- basename(env$OUTPUT_DIR)
    projectDir  <- env$OUTPUT_DIR
    getSamplePrefix <- function(sample) file.path(env$TASK_DIR, sample)
} else {
    projectName <- env$DATA_NAME
    projectDir  <- env$TASK_DIR
    getSamplePrefix <- function(sample) file.path(env$TASK_DIR, sample, sample)
}
#=====================================================================================

#=====================================================================================
# load the metadata and called SVs across all samples
#-------------------------------------------------------------------------------------
message("loading sample metadata")
inFile <- paste(env$FIND_PREFIX, "metadata", "yml", sep=".")
metadata <- read_yaml(inFile)
metadata <- lapply(metadata, function(x) strsplit(as.character(x), "\\s+")[[1]])
#=====================================================================================

#=====================================================================================
# load fixed-width bins data
#-------------------------------------------------------------------------------------
message("loading bin metadata")
binsDir <- file.path(env$GENOMES_DIR, "bins", env$GENOME, "fixed_width_bins")
binsFile <- paste(env$GENOME, "bins.size_65536.k_100.e_1.bed.gz", sep = ".")
binsFile <- file.path(binsDir, binsFile)
bins <- fread(
    binsFile,
    sep = "\t",
    header = FALSE,
    colClasses = c(
        "character", 
        "integer", 
        "integer", 
        "integer", 
        "numeric", 
        "character", 
        "integer", 
        "integer", 
        "integer", 
        "numeric", 
        "numeric"
    ),
    col.names = c(
        "chrom", 
        "start", 
        "end", 
        "cumIndex", 
        "gc", 
        "strand", 
        "excluded", 
        "gap", 
        "bad", 
        "umap", 
        "genmap"
    )
)
bins <- bins[, .SD, .SDcols = c(
    "chrom",
    "start",
    "gc",
    "excluded",
    "genmap"
)]
#=====================================================================================

#=====================================================================================
# create the composite bin coverage file
#-------------------------------------------------------------------------------------
message("merging sample coverage")
for(sample in metadata$SAMPLES){
    message(paste("   ", sample))
    coverageFile <- paste(getSamplePrefix(sample), env$GENOME, "coverage.index.gz", sep = ".")
    y <- fread(
        coverageFile,
        sep = "\t",
        header = TRUE,
        colClasses = c(
            "integer", # chromIndex, 1-referenced
            "integer", # chunkIndex, 0-referenced
            "integer", # cumNBreaks
            "numeric"  # coverage
        )
    )
    maxChunkIndices <- y[, .(maxChunkIndex = max(chunkIndex)), by = chromIndex]
    maxChunkIndices <- maxChunkIndices[order(chromIndex), c(0, cumsum(maxChunkIndex))]
    y[, ":="(
        chrom = metadata$CHROMS[chromIndex],
        start = chunkIndex * 65536 # chrom coordinate, 0-referenced (like BED)
    )]
    bins <- merge(
        bins,
        y[, .SD, .SDcols = c(
            "chrom",
            "start",
            "coverage"
        )],
        by = c("chrom", "start"),
        all.x = TRUE
    )
    names(bins)[ncol(bins)] <- sample
}
outFile <- paste(env$COVERAGE_PREFIX, 'rds', sep = ".")
saveRDS(bins, file = outFile)
#=====================================================================================
