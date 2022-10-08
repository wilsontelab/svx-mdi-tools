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
classDir <- file.path(env$MODULES_DIR, 'classes/R/nbinomCountsGC2')
sourceScripts(classDir, c('nbinomCountsGC2_class', 'nbinomCountsGC2_methods'))
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
#=====================================================================================

#=====================================================================================
# use the median percentile for each window to establish its corrected model
#-------------------------------------------------------------------------------------
expandWindows <- function(cell_id, data, key){ # go from window values to bin values, repeating window value over all bins
    cell <- cells[[cell_id]] # NOT data!
    w <- cell$window_size
    rollingRanges <- rollingRanges[[paste("w", w, sep = "_")]]
    chroms <- rollingRanges[reference_window == TRUE, chrom]
    unlist(rollingRanges[,
        sapply(data[[cell_id]][[key]][chroms == chrom[1]], function(x) rep(x, w))[1:.N],
        by = chrom
    ][, 2])
}
collapseWindows <- function(window_size, values){ # reverse the process of expandWindows
    chroms <- rowRanges[, chrom]
    unlist(rowRanges[,
        {
            i <- seq(1 + (window_size - 1) / 2, .N, window_size)
            values[chroms == chrom[1]][i]
        },
        by = chrom
    ][, 2])
}

# calculate the median percentile for each from among the bins and samples declared to be trustworthy
message("calculating bin median percentiles")
percentiles_w <- as.data.table(
    mclapply(constants$good_cell_ids, expandWindows, cells, "percentile", mc.cores = env$N_CPU) 
)
medianPercentile_w <- apply(percentiles_w, 1, median, na.rm = TRUE)
medianPercentile_wr <- mclapply(window_sizes, collapseWindows, medianPercentile_w, mc.cores = env$N_CPU) # 
names(medianPercentile_wr) <- paste("w", window_sizes, sep = "_")

# calculate the error correction per cell-window and re-run a modified HMM
message("recalculating good cells with batch effect correction")
goodCells <- mclapply(constants$good_cell_ids, function(cell_id){
    cell <- cells[[cell_id]]    

    # calculate the revise read expectation, i.e., reads per allele
    size <- cell$theta
    mu   <- cell$ER_gc / env$PLOIDY
    ww <- paste("w", cell$window_size, sep = "_")
    medianPercentile <- medianPercentile_wr[[ww]]
    excess <- qnbinom(medianPercentile, size = size, mu = mu) -
              qnbinom(0.5,              size = size, mu = mu) # the systematic bias across all cells       
    rpa <- mu + excess # thus, expect more reads for bins with consistently high nbinom percentile

    # run the final HMM
    rollingRanges <- rollingRanges[[ww]]
    hmm <- viterbi(
        cell$fit,
        binCounts = cell$NR_map_wr, 
        fractionGC = cell$gc_wr,
        percentile = medianPercentile,
        transProb = env$TRANSITION_PROBABILITY,
        chroms = rollingRanges[reference_window == TRUE, chrom],
        asRle = FALSE
    )    

    # return copy number profiles for plotting and clustering
    data.table(
        cn = cell$NR_map_wr / rpa,        
        hmm = hmm$cn
    )
}, mc.cores = env$N_CPU)
names(goodCells) <- constants$good_cell_ids
#=====================================================================================

#=====================================================================================
# calculate further parameters based on the fit across all cells
#-------------------------------------------------------------------------------------

# calculate the modal copy number for each bin
message("calculating bin modal copy number")
hmm_w <- as.data.table(
    mclapply(constants$good_cell_ids, expandWindows, goodCells, "hmm", mc.cores = env$N_CPU)
)
modalCN_w <- apply(hmm_w, 1, median, na.rm = TRUE)
modalCN_wr <- mclapply(window_sizes, collapseWindows, modalCN_w, mc.cores = env$N_CPU)
names(modalCN_wr) <- paste("w", window_sizes, sep = "_")

# infer sample sex from bin modal copy number
message("inferring sample sex")
X <- as.integer(round(median(modalCN_w[rowRanges$chrom == "chrX"], na.rm = TRUE), 0))
Y <- as.integer(round(median(modalCN_w[rowRanges$chrom == "chrY"], na.rm = TRUE), 0))
constants$sex <- paste0( 
    rep("X", if(is.na(X)) 0 else X), 
    rep("Y", if(is.na(Y)) 0 else Y),
    collapse = ""
)
message(paste(" ", constants$sex))

# calculate the copy number change per cell
message("calculating cell copy number change and finding CNVs")
for(cell_id in constants$good_cell_ids){
    cell <- goodCells[[cell_id]] 
    ww <- paste("w", cell$window_size, sep = "_")  
    goodCells[[cell_id]][, ":="(
        cnc = cell$cn  - modalCN_wr[[ww]],
        cnv = cell$hmm - modalCN_wr[[ww]]
    )]  
}
#=====================================================================================

#=====================================================================================
# save the final cell-level output file
# TODO: create a SummarizedExperiment object, perhaps as hdf5?
#-------------------------------------------------------------------------------------
message("saving final data object")
saveRDS(list(
    env = env,
    constants = constants,
    metadata = metadata,
    rowRanges = cbind( # identical for genome + maxMaxWindow
        rowRanges, 
        as.data.table(lapply(rollingRanges, function(x) x$reference_window))
    ), 
    colData = colData,
    goodCells = goodCells,  
    badCells = cells[constants$bad_cell_ids]
), file  = env$NORMALIZE_FILE)
#=====================================================================================
