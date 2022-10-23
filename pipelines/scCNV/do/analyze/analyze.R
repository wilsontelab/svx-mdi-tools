# analyze cells one at a time to set windows and fit to GC bias

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
message("initializing")
suppressPackageStartupMessages({
    library(data.table)
    library(zoo)
    library(parallel)
})
#-------------------------------------------------------------------------------------
# load the extract step file; do first to not overwrite any code below
load(Sys.getenv('EXTRACT_FILE'))
#-------------------------------------------------------------------------------------
# parse environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'PLOTS_DIR',
        'PLOT_PREFIX',
        'ANALYZE_FILE'
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
dir.create(env$PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
plotFilePrefix <- paste(env$PLOT_PREFIX, "gc_bias", sep = ".")
unlink(paste(plotFilePrefix, "*", sep = "."))
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
options(warn=2)
#-------------------------------------------------------------------------------------
# set some constants
ploidyFactor <- (env$PLOIDY + 0.5) / env$PLOIDY
minBinCount <- (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
#=====================================================================================

#=====================================================================================
# characterize the individual cells
#-------------------------------------------------------------------------------------
message('characterizing individual cells')
windowBins <- rowRanges$autosome & rowRanges$mappability >= env$MIN_MAPPABILITY
fitCell <- function(cell_id, minBinCount, pass){ # called twice: 1) learn about cell 2) optimize window size based on overdispersion

    # set the per-cell window size as the number of bins needed to obtain a median raw count >=minBinCount
    NR_raw_b <- raw_counts[[cell_id]][windowBins] # here, only use mappable autosome bins
    NR_med_b <- median(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_med_b) || NR_med_b == 0) return(NA)
    window_size <- ceiling(minBinCount / NR_med_b)
    window_size <- window_size + !(window_size %% 2) # make odd to ensure window centering

    # permanently remove cells with insufficient coverage, as determined by MAX_WINDOW_BINS
    if(window_size > env$MAX_WINDOW_BINS) return(NA)

    # calculate window sums and correct for mappability
    NR_raw_b <- raw_counts[[cell_id]] 
    rollingRanges <- rollingRanges[[paste("w", window_size, sep = "_")]]
    mappability <- rollingRanges[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_map_w <- unname(unlist(sapply(constants$chrom, function(chrom){
        chromIs  <- rowRanges$chrom == chrom
        NR_raw_w <- rollsum(NR_raw_b[chromIs], window_size, na.pad = TRUE, align = "center")
        NR_raw_w / mappability[chromIs]
    })))

    # process required bin and cell information
    reference_windows <- rollingRanges$reference_window # logical, never NA
    gc_w   <- rollingRanges$gc_fraction
    gc_wr  <- gc_w[reference_windows]
    gc_wra <- gc_w[reference_windows & rowRanges$autosome]
    NR_map_wr  <- NR_map_w[reference_windows] # THESE CAN BE NA DUE TO MABBABILITY
    NR_map_wra <- NR_map_w[reference_windows & rowRanges$autosome]
    chroms_wr <- rollingRanges$chrom[reference_windows]
    autosomes <- rollingRanges$chrom[reference_windows & rowRanges$autosome]

    # perform an initial GC bias fit that assumes the same CN for all autosomes
    fit <- new_nbinomCountsGC2(NR_map_wra, gc_wra, binCN = env$PLOIDY)

    # use the initial fit to solve an initial CN estimate for all autosomes
    hmm <- viterbi(fit, NR_map_wra, gc_wra, asRle = FALSE,
                   chroms = autosomes, transProb = env$TRANSITION_PROBABILITY)
    cn_estimate <- ifelse(hmm$cn == hmm$maxCN || hmm$cn == 0, NA, hmm$cn) # bins at maxCN are unreliable as they might be >maxCN

    # revise to a final GC bias fit using the initial copy number estimates
    fit <- new_nbinomCountsGC2(NR_map_wra, gc_wra, binCN = cn_estimate)

    # solve a final CN estimate for all chromosomes (not just autosomes)
    hmm <- viterbi(fit, NR_map_wr, gc_wr, asRle = FALSE, 
                   chroms = chroms_wr, transProb = env$TRANSITION_PROBABILITY)
    ER_gc <- predict(fit, gc_wr, type = 'adjustedPeak') * env$PLOIDY # use peak for visualization, unless it is unreliable
    if(pass == 1){
        cn <- NR_map_wr / ER_gc * env$PLOIDY
        percentile <- NULL
        theta <- NULL
    } else {
        cn  <- NULL # keep object size down
        percentile <- cumprob(fit, NR_map_wr, gc_wr, binCN = hmm$cn) 
        theta <- predict(fit, gc_wr, type = 'theta')
        hmm <- NULL
    }

    # return our results
    list(
        # common values
            rejected = TRUE, # caller must set to FALSE if keeping the cell
            pass = pass,
            window_size = window_size,
        # used by cells mclapply
            gc_wr = gc_wr,
            NR_map_wr = NR_map_wr,
            hmm = hmm, # in final pass, includes same windows as NR_map 
            cn = cn,
        # used by normalize
            fit = fit, 
            ER_gc = ER_gc,                       
            percentile = percentile,
            theta = theta
    )
}
cells <- mclapply(cell_ids, function(cell_id){
    cellPlotPrefix <- paste(plotFilePrefix, cell_id, sep = ".")
    
    # first pass fit at the sensitivity expected for Poisson without over/under-dispersion
    x1 <- fitCell(cell_id, minBinCount, pass = 1)
    if(!is.list(x1)) return(list(
        rejected = TRUE,
        pass = 1,
        window_size = NA # the very worst cells, wholly insufficent data
    ))
    plot(x1$fit, x1, env$PLOIDY, cellPlotPrefix)

    # second pass fit at a sensitivity adjusted for the specific cell's dispersion
    dcn <- x1$cn[x1$hmm$cn == env$PLOIDY] - env$PLOIDY
    scaleFactor <- env$N_SD_HALFCN * sd(dcn, na.rm = TRUE) / 0.5
    x <- fitCell(cell_id, minBinCount * scaleFactor, pass = 2)
    if(!is.list(x)) return({
        x1$cn  <- NULL
        x1$hmm <- NULL
        x1
    })
    x$rejected <- FALSE    
    plot(x$fit, x, env$PLOIDY, cellPlotPrefix)  
    x
}, mc.cores = env$N_CPU)
names(cells) <- cell_ids

# assemble and organize the per-cell data
message('recording cell metadata')
colData[, window_size := 0L] # needed to prevent errors if first cell has window_size == NA
colData[, ':='(
    rejected    = cells[[cell_id]]$rejected,
    pass        = cells[[cell_id]]$pass,
    window_size = cells[[cell_id]]$window_size
), by = cell_id]

#=====================================================================================
# save the interim output file
#-------------------------------------------------------------------------------------
rm(raw_counts)
save.image(file = env$ANALYZE_FILE)


saveRDS(list(
    plotFilePrefix = plotFilePrefix,
    colData = colData,
    cells = cells
), file  = env$ANALYZE_RDS)

#=====================================================================================
