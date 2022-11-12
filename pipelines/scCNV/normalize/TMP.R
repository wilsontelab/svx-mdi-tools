#=====================================================================================
# extract and reformat data from step 'bin', or 10x Cell Ranger DNA, as input
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
        'EXTRACT_FILE'        
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
sourceScripts(scCnvSharedDir, c('loadHdf5', 'parseGenomeBins'))
#=====================================================================================



# determine the number of reads per bin required for a cell's data to support a robust HMM
# with default options, places ~96% of bins within modal_CN +/- 0.5
getMinBinCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
    minBinCount <- (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# set the per-cell window size as the number of bins needed to obtain a median raw count >=minBinCount
# permanently remove cells with insufficient coverage, as determined by MAX_WINDOW_BINS
getCellWindowSize <- function(cell_id, minBinCount, bins = NULL){
    return(1)
}

# fit a cell's GC bias using the negative binomial distribution
# called twice: 1) learn about cell 2) optimize window size based on overdispersion
fitCell_ <- function(cell_id, stage = "extract", pass = 1, 
                     scaleFactor = NULL, prevFit = NULL){ 

    # calculate window sums and correct for mappability
    modal_CN <- colData[cell_id, modal_CN]
    minBinCount <- if(is.null(scaleFactor)) getMinBinCount(modal_CN) else scaleFactor * prevFit$minBinCount
    window_size <- getCellWindowSize(cell_id, minBinCount)
    if(is.na(window_size)) return(prevFit)    
    if(!is.null(prevFit) && window_size == prevFit$window_size) return(prevFit)
    windows <- windows[[paste("w", window_size, sep = "_")]]
    mappability <- windows[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_map_w <- unname(unlist(sapply(constants$chrom, function(chrom){
        collapseVector(
            raw_counts[[cell_id]][rowRanges$chrom == chrom], 
            window_size
        ) / mappability[windows$chrom == chrom]
    })))

    # process required bin and cell information
    gc_w  <- windows$gc_fraction
    gc_wa <- gc_w[windows$autosome]
    NR_map_wa <- NR_map_w[windows$autosome]
    autosomes <- windows$chrom[windows$autosome]

    # perform an initial GC bias fit that assumes the same CN for all autosomes
    fit <- new_nbinomCountsGC2(NR_map_wa, gc_wa, binCN = modal_CN)

    # use the initial fit to solve an initial CN estimate for all autosomes
    hmm <- viterbi(fit, NR_map_wa, gc_wa, asRle = FALSE,
                   chroms = autosomes, transProb = env$TRANSITION_PROBABILITY)
    cn_estimate <- ifelse(hmm$cn == hmm$maxCN | hmm$cn == 0, NA, hmm$cn) # bins at maxCN are unreliable as they might be >maxCN

    # revise to a final GC bias fit using the initial copy number estimates
    fit <- new_nbinomCountsGC2(NR_map_wa, gc_wa, binCN = cn_estimate)

    # solve a final CN estimate for all chromosomes (not just autosomes)
    ER_gc <- predict(fit, gc_w, type = 'adjustedPeak') * modal_CN # use peak for visualization, unless it is unreliable   
    cn <- NR_map_w / ER_gc * modal_CN 
    hmm <- viterbi(fit, NR_map_w, gc_w, asRle = FALSE, 
                   chroms = windows$chrom, transProb = env$TRANSITION_PROBABILITY)

    # return our results
    list(
        rejected = TRUE, # caller must set to FALSE if keeping the cell
        stage = stage,
        pass = pass,
        modal_CN = modal_CN,
        minBinCount = minBinCount,
        window_size = window_size,
        NR_map_w = NR_map_w,            
        gc_w = gc_w,
        fit = fit, 
        ER_gc = ER_gc,  
        cn = cn,              
        hmm = hmm$cn,
        percentile = cumprob(fit, NR_map_w, gc_w, binCN = hmm$cn),         
        theta = predict(fit, gc_w, type = 'theta')
    )
}
fitCell <- function(cell_id, stage = "extract"){ 

    # the very worst cells, wholly insufficent data
    x0 <- list(
        rejected = TRUE,
        stage = stage, 
        pass = 0,
        window_size = NA 
    )

    # first pass fit at a sensitivity adjusted for a specific cell's read depth
    x1 <- fitCell_(cell_id, stage, pass = 1)
    if(!is.list(x1)) return(x0)
    x1$rejected <- FALSE
    x1


    # # second pass fit at a sensitivity adjusted for the specific cell's read depth and dispersion
    # dcn <- x1$cn[x1$hmm == x1$modal_CN] - x1$modal_CN
    # scaleFactor <- max(1, env$N_SD_HALFCN * sd(dcn, na.rm = TRUE) / 0.5, na.rm = TRUE)
    # x2 <- if(scaleFactor == 1) x1
    #       else fitCell_(cell_id, stage, pass = 2, scaleFactor, x1)
    # dcn <- x2$cn[x2$hmm == x2$modal_CN] - x2$modal_CN
    # sdDcn <- sd(dcn, na.rm = TRUE)
    # x2$rejected <- is.na(sdDcn) || sdDcn * env$N_SD_HALFCN > 1 # deliberately 2-fold more permissive than the stated target
    # x2
}

# make a composite plot of cell model for QC purposes
plotCellByWindow <- function(d, ylab, ylim = NULL){
    plot(1:length(d), d, bty = "n",
        pch = 19, cex = 0.4, 
        col = rgb(0, 0, 0, 0.2),
        xaxt = "n", xlab = NULL, ylab = ylab, ylim = ylim)
}
plotCellQC <- function(cell_id, cell){
    png(
        filename = paste(env$PLOT_PREFIX, cell_id, "qc", "png", sep = "."),
        width  = 1.5 * 6, 
        height = 1.35, 
        units = "in", 
        pointsize = 7,
        bg = "white",
        res = 96, # i.e., optimized for on-screen display
        type = "cairo"
    )
    layout(matrix(c(c(1,1), rep(c(2,3), 5)), nrow = 2, ncol = 6))

    # plot NR_map_w vs. gw_w
    par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
    plot(cell$fit, cell$gc_w, cell$NR_map_w, cell$modal_CN, cell$rejected)

    # plot NR_map_w vs. window index, i.e., pre-normalization
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plotCellByWindow(cell$NR_map_w, "# Reads")

    # plot CN vs. window index, i.e., post-normalization
    plotCellByWindow(cell$cn, "CN", ylim = c(0, 4))
    abline(h = 0:4, col = "grey")
    if(!is.null(cell$hmm)) lines(1:length(cell$hmm), cell$hmm, col = "red")

    dev.off()
}



#=====================================================================================
# characterize the individual cells
# first pass fit at the bin resolution expected for Poisson without over/under-dispersion
#-------------------------------------------------------------------------------------
message('characterizing individual cells')
colData[, modal_CN := env$PLOIDY]

# # 27026 APC
# ALLOWED <- as.character(c(43,34,38,17,41,27,16,29,19,32,4,35,26,33,6))
# 26898_KRAS
ALLOWED <- as.character(c(13,3,17,7,0,5,10,1,12,6,11))

cells <- mclapply(cell_ids, function(cell_id){
    if(!(cell_id %in% ALLOWED)) return(list(
        rejected = TRUE,
        stage = "extract", 
        pass = 0,
        window_size = NA 
    ))
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
    cells = cells
), file  = env$EXTRACT_FILE)
#=====================================================================================
