# determine the least number of reads per bin for a cell's data to support a robust HMM
# with default options, places ~96% of bins within modal_CN +/- 0.5
getMinBinCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
    minBinCount <- (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# set the per-cell window size as the number of bins needed to obtain a mean raw count >=minBinCount
# permanently remove cells with insufficient coverage, as determined by MAX_WINDOW_BINS
getMinWindowPower <- function(cell_id, minBinCount, bins = NULL){
    if(is.null(bins)) bins <- rowRanges$autosome & rowRanges$mappability >= env$MIN_MAPPABILITY
    NR_raw_b <- raw_counts[[cell_id]][bins]
    NR_avg_b <- mean(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_avg_b) || NR_avg_b == 0) return(NA)
    window_size <- ceiling(minBinCount / NR_avg_b)
    if(is.na(window_size) || window_size > 2 ** env$MAX_WINDOW_BINS) return(NA)
    window_sizes <- 2 ** (0:env$MAX_WINDOW_BINS)
    window_size <- window_sizes[min(which(window_sizes >= window_size))] # thus, allow 2**(0:MWB) bins per window
    log2(window_size)
}

# fit a cell's GC bias using the negative binomial distribution
# called as many times as needed at increasing window sizes to account for a cell's overdispersion
fitCell_ <- function(cell_id, modal_CN, minBinCount, windowPower, pass, stage = "extract"){ 

    # calculate window sums and correct for mappability
    window_size <- 2 ** windowPower
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
    x <- list(
        rejected = TRUE,
        stage = stage, 
        pass = 0,
        window_size = NA 
    )

    # get the smallest reasonable window size for this cell based on its read depth, assuming Poisson
    modal_CN <- colData[cell_id, modal_CN]
    minBinCount <- getMinBinCount(modal_CN)
    minWindowPower <- getMinWindowPower(cell_id, minBinCount)
    if(is.na(minWindowPower)) return(x)

    # step up to larger and larger windows as needed to account for overdispersion
    # expanding window sizes by a factor of two leads to predictable cell-to-cell bin relationships
    pass <- 0
    windowPower <- minWindowPower - 1 
    while(x$rejected && windowPower < env$MAX_WINDOW_BINS){
        pass <- pass + 1
        windowPower <- windowPower + 1
        x <- fitCell_(cell_id, modal_CN, minBinCount, windowPower, pass, stage)
        dcn <- x$cn[x$hmm == x$modal_CN] - x$modal_CN
        sdDcn <- sd(dcn, na.rm = TRUE)
        x$rejected <- is.na(sdDcn) || sdDcn * env$N_SD_HALFCN > 1 # deliberately 2-fold more permissive than the stated target
    }

    # return our result
    x
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
