
# create a sample manifest for use by the svx family of apps

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
        'PIPELINE_NAME',
        'OUTPUT_DIR',
        'TASK_DIR',
        'DATA_NAME',
        'FIND_PREFIX',
        'MANIFEST_PREFIX',
        'FIND_MODE'
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
} else {
    projectName <- env$DATA_NAME
    projectDir  <- env$TASK_DIR
}
#=====================================================================================

#=====================================================================================
# load the metadata and called SVs across all samples
#-------------------------------------------------------------------------------------
message("loading sample metadata")
inFile <- paste(env$FIND_PREFIX, "metadata", "yml", sep=".")
metadata <- read_yaml(inFile)
metadata <- lapply(metadata, function(x) strsplit(x, "\\s+")[[1]])

message("loading SV summary table")
inFile <- paste(env$FIND_PREFIX, 'structural_variants', 'gz', sep = ".")
svs <- fread(inFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#=====================================================================================

#=====================================================================================
# create the sample manifest table
#-------------------------------------------------------------------------------------
message("writing sample manifest")
getSvxStatValue <- function(action, scriptId, statName){
    mapply(getStatValue, 
           projectDir, metadata$SAMPLES, 
           action, scriptId, statName)
}
manifest <- data.table(
    Project             = projectName, 
    Sample_ID           = metadata$SAMPLES, 
    Description         = metadata$SAMPLES,
    maxTLen             = as.integer(metadata$MAX_TLENS),
    nReadPairs          = getSvxStatValue('align', 'prepare_fastq', 'nInputPairs'),
    nSourceMolecules    = getSvxStatValue('extract', 'extract_nodes', 'nReadPairs'),    
    onTargetCoverage    = getSvxStatValue('collate', 'coverage', 'foldCoverage[TT]'),
    offTargetCoverage   = getSvxStatValue('collate', 'coverage', 'foldCoverage[--]'),
    enrichment          = getSvxStatValue('collate', 'coverage', 'enrichment'),
    nSvs                = sapply(metadata$SAMPLES, function(x) sum(svs[[x]] > 0))
)
manifest[, Yield := nReadPairs]
outFile <- paste(env$MANIFEST_PREFIX, 'sample_manifest', 'csv', sep = ".")
fwrite(
    manifest, 
    file = outFile, 
    quote = FALSE, 
    sep = ",",    
    row.names = FALSE,   
    col.names = TRUE
)
#=====================================================================================
