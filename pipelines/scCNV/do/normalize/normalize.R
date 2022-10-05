# perform cross-cell normalization, recall single-cell CNVs

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
message("initializing")
suppressPackageStartupMessages({
    # library(rhdf5)
    library(data.table)
    # library(SummarizedExperiment)
    # library(GenomicRanges)
    library(zoo)
    library(parallel)
})
#-------------------------------------------------------------------------------------
# load the extract step file; do first to not overwrite any code below
load(Sys.getenv('ANALYZE_FILE'))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'NORMALIZE_FILE'
    ),
    integer = c(
        "PLOIDY",
        "N_CPU"
    ),
    double = c(
        "TRANSITION_PROBABILITY"
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
classDir <- file.path(env$MODULES_DIR, 'classes/R/nbinomCountsGC')
sourceScripts(classDir, c('nbinomCountsGC_class', 'nbinomCountsGC_methods'))
classDir <- file.path(env$MODULES_DIR, 'classes/R/hmmEPTable')
sourceScripts(classDir, c('hmmEPTable_class', 'hmmEPTable_methods'))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# set some "constants", adding to values placed by Cell Ranger
colData <- rbind(colData[rejected == FALSE], colData[rejected == TRUE])
constants$good_cell_ids  <- colData[rejected == FALSE, cell_id]
constants$bad_cell_ids   <- colData[rejected == TRUE,  cell_id]
constants$num_good_cells <- length(constants$good_cell_ids)
constants$num_bad_cells  <- length(constants$bad_cell_ids)
constants$num_good_bins  <- nrow(rowRanges)
constants$num_bad_bins   <- constants$genome_bins - constants$num_good_bins
# TODO: determine sample sex? probably do this below using percentiles
#=====================================================================================

#=====================================================================================
# use the median percentile for each window to establish its corrected model
#-------------------------------------------------------------------------------------

# calculate the median percentile for each from among the bins and samples declared to be trustworthy
message("calculating median percentile")
percentiles <- as.data.table(lapply(constants$good_cell_ids, function(cell_id){
    cellFits[[cell_id]]$percentile
}))
medianPercentile <- apply(percentiles, 1, median, na.rm = TRUE)

# calculate the error correction per cell-window and re-run a modified HMM
message("recalculating cells with batch effect correction")
normalized <- mclapply(constants$good_cell_ids, function(cell_id){
    cell <- cellFits[[cell_id]]    

    # calculate the revise read expectation, i.e., reads per allele
    size <- cell$fit$theta
    mu <- cell$ER_gc / env$PLOIDY
    excess <- qnbinom(medianPercentile, size = size, mu = mu) -
              qnbinom(0.5,              size = size, mu = mu) # the systematic bias across all cells       
    rpa <- mu + excess # thus, expect more reads for bins with consistently high nbinom percentile

    # run the final HMM
    rollingRanges <- rollingRanges[[paste("w", cell$window_size, sep = "_")]]
    hmm <- viterbi(
        cell$fit,
        binCounts = cell$NR_map_w, 
        fractionGC = rollingRanges$gc_fraction,
        percentile = medianPercentile,
        transProb = env$TRANSITION_PROBABILITY,
        chroms = rollingRanges$chrom,
        asRle = FALSE
    )    

    # return copy number profiles for plotting and clustering
    list(
        hmm = hmm$cn,
        cn = cell$NR_map_w / rpa
    )
}, mc.cores = env$N_CPU)
names(normalized) <- constants$good_cell_ids
#=====================================================================================

#=====================================================================================
# calculate further parameters based on the fit across all cells and assemble the final layers
#-------------------------------------------------------------------------------------
message("packaging final data objects")
buildCellLayer <- function(key) as.data.table(lapply(constants$good_cell_ids, function(cell_id){
    normalized[[cell_id]][[key]]
}))
cn  <- buildCellLayer("cn")
hmm <- buildCellLayer("hmm")
bad <- lapply(constants$bad_cell_ids, function(cell_id) cellFits[[cell_id]]$cn)
setnames(cn,  constants$good_cell_ids)
setnames(hmm, constants$good_cell_ids)
names(bad) <- constants$bad_cell_ids
rowRanges[, cn := apply(hmm[, .SD, .SDcols = constants$good_cell_ids], 1, median, na.rm = TRUE)]
rollingRanges <- as.data.table(lapply(rollingRanges, function(x) x$reference_window))
#=====================================================================================

#=====================================================================================
# save the final cell-level output file
# TODO: create a SummarizedExperiment object, perhaps as hdf5?
#-------------------------------------------------------------------------------------
message("saving final data objects")
pipelineEnv <- env
save(
    pipelineEnv,
    constants, 
    metadata,
    rowRanges,
    rollingRanges,
    colData,
    cn,
    hmm,
    bad,
    file = env$NORMALIZE_FILE
)
#=====================================================================================
