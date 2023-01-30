#=====================================================================================
# normalize and analyze individual cells
#-------------------------------------------------------------------------------------
# throughout this script and pipeline, the following naming conventions are used
#-------------------------------------------------------------------------------------
#   variable names:
#     NR  = Number of Reads, i.e., a read count (can be fractional due to fragment splitting across bins)
#     ER  = Expected Reads, i.e., the read count that corresponds to a specific copy number and/or state
#     RPA = Reads Per Allele, i.e., the read count expected at CN == 1, i.e., ER_ploidy / ploidy
#     CN  = Copy Number, PRIOR to any replication of an allele, can be fractional estimates as NR / ER
#     HMM = Hidden Markov Model, i.e., quantal CN values determined from a model of NR + ER
#     NA_ = Number of Alleles, i.e., how many DNA copies are present considering BOTH CN and replication (underscore since NA is reserved)
#     NAR = Number of Alleles Replicated, thus, could be c(0,1,2) for CN==2, to yield NA_=c(2,3,4)
#     FAR = Fraction of Alleles Replicated, i.e., sum(NAR) / sum(CN), synonymous with fractionS
#-------------------------------------------------------------------------------------
#   variable "subscripts", e.g., XX_xx:
#     g = a value or list applied to an entire Genome
#     c = a value or list applied to an entire Chromosome
#     b = a value applied to 20kb Bins
#     w = a value applied to Windows, i.e., sequential groups of input bins
#     m = counts that have been subjected to a Mappability correction
#     f = bins or windows that have been Filtered, to remove some bins during downstream analysis
#     l = bins or windows that have been subjected to Limits, e.g., to exlude outlier windows
#     d = bins or windws that have been Downsampled to improve fitting speed
#     h = counts that have been rescaled by HMM toward ploidy
#     r = counts that have been rescaled by NAR toward unreplicated
#     s = counts that have been rescaled based on chromosome shape biases (see below)
#-------------------------------------------------------------------------------------
# additional considerations:
#   - "ploidy" is the expected autosome copy number prior to replication, either 1 or 2
#   - "modal" and "peak" are largely interchangeable and refer to the mode of a set of NR_w values
#   - "windowSize" is the number of 20kb bins that comprise a cell's windows, can vary between cells
#   - "windowPower" = log2(windowSize), i.e., allowed window sizes increment by doubling
#   - "GC bias" is the dependence of NR values on window GC base content
#   - "squashed" refers to a GC bias correction that does not consider replication, 
#                 i.e., that masks replication-based differences in NR as a result of GC bias correction
#   = "shape" is an artifact where a chromosome's counts are not flat but tilted, S-shaped, etc.
#                 shape effects can vary considerably with the amplification method
#   - "fractionS" is the total fraction of alleles across the genome that have been replicated (FAR), 
#                 e.g., at fractionS==0.5, half of the genome's bases (not windows!) have been replicated
#=====================================================================================

#=====================================================================================
# some general constants and working variables, others also found below in context
#-------------------------------------------------------------------------------------
workingStep <- "initializing"
maxChromCN <- 5 # thus, CN==5 is not trustworthy, the true value could be higher
chromCNs <- 0:maxChromCN
maxSequentialCN <- 5
maxCompositeCN <- 4
repTransProb <- 1e-1
#=====================================================================================

#=====================================================================================
# generic numeric and distribution support functions
#-------------------------------------------------------------------------------------
peakValue <- function(x){ # find the peak, i.e., mode of a set of values
    x <- x[!is.na(x)]
    if(length(x) == 0) return(NA)
    d <- density(x)
    d$x[which.max(d$y)]
}
excludeOutliers <- function(v, low = 0.025, high = 0.975, min = -Inf){ # could be more sophisticated, but these are heuristics
    v <- v[!is.na(v)]
    q <- quantile(v, c(low, high))
    v[v >= q[1] & v <= q[2] & v >= min]
}
#=====================================================================================

#=====================================================================================
# plotting functions, normally only used during development and debugging
#-------------------------------------------------------------------------------------
plotNumber <- 1
pointOpacity <- 0.15
defaultPointColor <- rgb(0, 0, 0, pointOpacity)
plotsDir <- file.path(env$PLOTS_DIR, 'scCNV')
if(!dir.exists(plotsDir)) dir.create(plotsDir)
####################
# plotsDir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV/plots"
saveCellPlot <- function(cell, name, fn, width = 2, height = 2, increment = TRUE){ # save a plot for future assembly and app
    if(!env$VERBOSE_PLOTS) return(NULL)
    plotDir <- file.path(plotsDir, cell$cell_id)
    if(!dir.exists(plotDir)) dir.create(plotDir)
    filename <- paste(cell$cell_id, sprintf("%02d", plotNumber), shapeKey, name, "png", sep = ".")
    pngFile <- file.path(plotDir, filename)
    png(pngFile, width = width, height = height, units = "in", pointsize = 7, res = 300, type = "cairo")
    par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
    fn()
    dev.off()
    if(increment) plotNumber <<- plotNumber + 1
}
getRepColor <- function(cell, default = NULL){ # after replication fitting, color points by replication state
    if(is.null(cell$replicationModel) || is.null(cell$replicationModel[[shapeKey]])) {
        if(is.null(default)) defaultPointColor else default
    } else {
        cols <- c(                            # by NAR (same order, different meaning, from CN below)
            rgb(0, 0, 1, pointOpacity),       # 0 = blue
            rgb(0.1, 0.8, 0.1, pointOpacity), # 1 = subdued green
            rgb(1, 0, 0, pointOpacity),       # 2 = red
            rgb(1, 0.65, 0, pointOpacity),    # 3 = orange
            defaultPointColor                 # 4 = yelblack/grey to make it obviouslow
        )
        repKey <- if(is.null(cell$windows[[shapeKey]]$composite)) "sequential" else "composite"
        col <- cols[cell$windows[[shapeKey]][[repKey]]$NAR + 1]
        col[is.na(col)] <- defaultPointColor
        col
    }
}
getCnColor <- function(cell){
    cols <- c(                            # by CN (regardless of ploidy)
        defaultPointColor,                # 0 = black/grey (absence of color/copies...)
        rgb(0, 0, 1, pointOpacity),       # 1 = blue ("cool" colors are losses)
        rgb(0.1, 0.8, 0.1, pointOpacity), # 2 = subdued green ("good/go" for typical CN neutral)
        rgb(1, 0, 0, pointOpacity),       # 3 = red ("hot" colors are gains)
        rgb(1, 0.65, 0, pointOpacity),    # 4 = orange
        defaultPointColor                 # 5 = back to black/grey to make it obvious
    ) 
    repKey <- if(is.null(cell$windows[[shapeKey]]$composite)) "sequential" else "composite"
    col <- cols[cell$windows[[shapeKey]][[repKey]]$HMM + 1]
    col[is.na(col)] <- defaultPointColor
    col
}
minGciColor <- 30
maxGciColor <- 55
gcPalette <- colorRampPalette(c(
    rgb(1, 0.65, 0, 0.3),
    rgb(0, 0,    1, 0.3)
), alpha = TRUE)(maxGciColor - minGciColor + 1)
getGcColor <- function(cell, windowI){
    gcis <- windows[[cell$windowPower + 1]][windowI, round(gc_fraction * 100, 0)]
    sapply(gcis, function(gci) {
        if(is.na(gci)) NA
        else if(gci < minGciColor) gcPalette[minGciColor - minGciColor + 1]
        else if(gci > maxGciColor) gcPalette[maxGciColor - minGciColor + 1]
        else gcPalette[gci - minGciColor + 1] 
    }) 
}
getWindowCoordinates <- function(cell, windowI){
    windows <- windows[[cell$windowPower + 1]]
    x <- (1:nrow(windows))[windowI]
    list(
        x = x,
        xlim = c(0, windows[, max(i) + 1]),
        chromV = c(0, windows[, max(i), by = "chrom"][[2]]) + 0.5,
        chromLab = gsub("chr", "", windows[, unique(chrom)]),
        chromLabX = windows[, {
            minI <- min(i)
            minI + (max(i) - minI) / 2
        }, by = "chrom"][[2]]
    )
}
addChromLabels <- function(coord){
    mtext(coord$chromLab, side = 1, line = ((1 + 1:length(coord$chromLab)) %% 2) + 0.5, 
          at = coord$chromLabX, cex = 1)
    abline(v = coord$chromV, col = "grey")
}
plotGcBias <- function(cell, name, composite, gc_w, NR_wm, windowCN = NULL, gc_fit = NULL, windowI = TRUE){ # NR vs. fraction GC, sometimes with adjustments for CN and/or replication
    saveCellPlot(cell, paste("gc", name, sep = "_"), function(){
        if(is.null(gc_fit)) gc_fit <- getGcFit(cell, composite)
        plot(gc_fit, gc_w, NR_wm, cell$ploidy, col = getRepColor(cell, NA)[windowI], binCN = windowCN)
    })
}
plotP_gc_nar <- function(cell, dt, iteration){ # NR vs. fraction GC, sometimes with adjustments for CN and/or replication
    saveCellPlot(cell, paste("P_gc_nar", iteration, sep = "_"), function(){
        dt <- dt[order(GCI)]
        dt$GC <- dt$GCI / nGcSteps
        ymax <- max(dt[["0"]], dt$all)
        plot(NA, NA, xlim = c(0.3, 0.6), ylim = c(0, ymax),
                     xlab = "Fraction GC", ylab = "Frequency" )
        cols <- list(
            "0" = "blue",
            "1" = "purple",
            "2" = "red3",
            "all" = "black"
        )
        for(nar in c("all", "0", "1", "2")){
            y <- dt[[nar]]
            if(!is.null(y)) lines(dt$GC, y, col = cols[[nar]])
        }
        legend("topright", legend = names(cols), col = unlist(cols), lwd = 1)
    })
}
plotWindows_counts <- function(cell, name, NR_wm, windowI = TRUE,  # plot window counts by chromosome coordinates
                               col = defaultPointColor, shape = NULL, 
                               increment = TRUE){
    saveCellPlot(cell, paste("NR_wm", name, sep = "_"), function(){
        coord <- getWindowCoordinates(cell, windowI)
             if(is.character(col) && col[1] == "useRepColor") col <- getRepColor(cell)
        else if(is.character(col) && col[1] == "useCnColor")  col <- getCnColor(cell)
        else if(is.character(col) && col[1] == "useGcColor")  col <- getGcColor(cell, windowI)
        else if(is.character(col) && col[1] == "badCell")     col <- rgb(1, 0, 0, pointOpacity)
        ER_modal_NA <- cell$windows[[shapeKey]]$ER_modal_NA
        plot(NA, NA, xlim = coord$xlim, ylim = c(0, ER_modal_NA * 3),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "# of Reads")
        addChromLabels(coord)
        RPA <- cell$windows[[shapeKey]]$RPA
        if(!is.null(RPA)) abline(h = RPA * 0:10, col = "grey")
        points(coord$x, NR_wm, pch = 19, cex = 0.3, col = col) # caller is responsible for applying shapes
        if(!is.null(shape)) points(coord$x, ER_modal_NA * shape, pch = 16, cex = 0.5, col = "blue")
    }, width = 6, increment = increment)
}
plotWindows_cn <- function(cell, name, composite, gc_w, NR_wm, windowI = TRUE, 
                           hmm = NULL, col = "replication", gc_fit = NULL, shape = NULL){ # plot window CN by chromosome coordinates
    saveCellPlot(cell, paste("cn", name, sep = "_"), function(){
        if(is.null(gc_fit)) gc_fit <- getGcFit(cell, composite)
        RPA_peak <- predict(gc_fit, gc_w, type = 'adjustedPeak')
        cn <- NR_wm / RPA_peak
        coord <- getWindowCoordinates(cell, windowI)
        col <- if(col == "CN") getCnColor(cell) 
               else if (col == "replication") getRepColor(cell)
               else if (col == "GC") getGcColor(cell, windowI)
               else col
        plot(NA, NA, xlim = coord$xlim, ylim = c(0, maxSequentialCN + 1),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "Copy Number")
        addChromLabels(coord)
        abline(h = 0:(maxSequentialCN + 1), col = "grey")
        points(coord$x, cn, pch = 19, cex = 0.3, col = col[windowI])
        if(!is.null(hmm)) lines(coord$x, hmm, col = "red3")
        if(!is.null(shape)) points(coord$x, cell$ploidy * shape, pch = 16, cex = 0.5, col = "blue")
    }, width = 6)
}
plotWindows_FAR_isolated <- function(cell, name, composite = FALSE){ # plot window CN by chromosome coordinates
    saveCellPlot(cell, paste("FAR", name, sep = "_"), function(){
        coord <- getWindowCoordinates(cell, TRUE)
        col <- getRepColor(cell)
        plot(NA, NA, xlim = coord$xlim, ylim = c(-0.5, 2),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "Fraction Replicated")
        addChromLabels(coord)
        abline(h = 0:1, col = "grey")
        repKey <- getRepModelKey(composite)
        points(coord$x, cell$windows[[shapeKey]][[repKey]]$FAR, pch = 19, cex = 0.3, col = col)
    }, width = 6)
}
plotWindows_CN_isolated <- function(cell, name, composite = FALSE){ # plot window CN by chromosome coordinates
    saveCellPlot(cell, paste("CN", name, sep = "_"), function(){
        coord <- getWindowCoordinates(cell, TRUE)
        col <- getCnColor(cell)
        plot(NA, NA, xlim = coord$xlim, ylim = c(0, maxSequentialCN + 1),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "Copy Number")
        addChromLabels(coord)
        abline(h = 0:(maxSequentialCN + 1), col = "grey")
        repKey <- getRepModelKey(composite)
        points(coord$x, cell$windows[[shapeKey]][[repKey]]$CN, pch = 19, cex = 0.3, col = col)
        lines(coord$x, cell$windows[[shapeKey]][[repKey]]$HMM, col = "red3")
    }, width = 6)
}
plotCellFit <- function(cell, name, model, bestModelByType){ # distributions of read counts with best-fitted replication models
    saveCellPlot(cell, paste("repFit", name, sep = "_"), function(){
        NR_a <- round(cell$tmp$NR_wmshfl, 0) # here, a is for "actual"
        N <- length(NR_a)
        NR_x <- 1:max(NR_a, na.rm = TRUE)
        getModelDensity <- function(model){ with(model, {
            weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, fractionS) # see utilities/create_fractionS_table.R
            dnbinom(NR_x, mu = ER_unreplicated * 1,   size = theta) * weights[1] +   
            dnbinom(NR_x, mu = ER_unreplicated * 1.5, size = theta) * weights[2] + 
            dnbinom(NR_x, mu = ER_unreplicated * 2,   size = theta) * weights[3]         
        })}
        NR_ad <- density(NR_a)
        plot(NR_ad, typ="l", col = "black", lwd = 1.5, main = "", xlab = "# of Reads",
             xlim = c(0, max(NR_ad$x)), ylim = c(0, max(NR_ad$y) * 1.1))
        abline(v = cell$windows[[shapeKey]]$ER_modal_NA * c(1/2, 1, 2))
        bestModelByType[, {
            density <- getModelDensity(.SD)
            col <- modelTypeColors[[modelType]]
            lwd <- if(modelType == model$modelType) 1.5 else 0.75
            lines(NR_x, density, col = col, lwd = lwd)
            abline(v = ER_unreplicated * 1:2, lty = 2, col = col)
        }, by = "modelType"]
    }, increment = FALSE)
}
plotCountDensity <- function(cell){ # this function plots the initial fits of nrModelTypes to NR_wm for different windowPowers
    saveCellPlot(cell, paste("NR_dist", paste0("windowPower", cell$windowPower), sep = "_"), function() with(cell$windows[[shapeKey]], {
        cnMax <- 5
        NR_min <- max(RPA / 10, 2)
        NR_max <- RPA * cnMax
        plot(NA, NA, xlab = "# of Reads", ylab = "density",
             xlim = c(NR_min, NR_max), ylim = c(0, max(density$density) * 1.1))
        abline(v = RPA * 1:5, col = "grey")
        with(density, { lines(NR, density, col = "black", lwd = 1.5) })
        col <- if(allow) "blue" else "red3"
        d <- nrModelTypes[[nrModelType]]$optFn(unlist(cell$windows[[shapeKey]][1:3]), density$NR)
        lines(density$NR, d, col = col)
        legend("topright", legend = nrModelType, col = col, lwd = 1, cex = 0.75)
    }), increment = FALSE)  
}
plotBatchEffect <- function(cell){
    saveCellPlot(cell, "batchEffect", function() {
        repKey <- getRepModelKey(cell$cellIsReplicating)
        getShapeDist <- function(shapeKey){
            density(cell$windows[[shapeKey]][[repKey]]$FAR, na.rm = TRUE)
        }
        d_s <- getShapeDist(if(env$SHAPE_CORRECTION %in% c('cell', 'both')) "shaped" else "unshaped")
        d_b <- getShapeDist("batched")
        plot(NA, NA, xlab = "Number of Alleles", ylab = "density",
             xlim = c(-1, 2), ylim = c(0, max(d_s$y, d_b$y) * 1.1))
        abline(v = seq(0, 1, 0.5), col = "grey")
        lines(d_s, col = "black")
        lines(d_b, col = "blue")
    }, increment = FALSE)  
}
#=====================================================================================

#=====================================================================================
# establish the common parameters of a set of mixture models to be applied to NR_wm distributions of individual cells
# the concept is that NR is quantal by NA at each window, regardless of the reason (replication, CNV, etc.)
# and that usually one one or two NR_wm peaks dominate in a cell (or sometimes three, with both chromosome losses and gains)
# by allowing for these different peak patterns, the best mu and theta values of the negative binomial distribution can be used
# to set to the target width (in CN units) of the peak at the modal CN to optimize CNV calling resolution and visualization
#-------------------------------------------------------------------------------------
nrModelTypes <- list(
    singlePeak = list( # all windows are sampled from a single distribution at NR = NR_ploidy
        par   = c(),
        lower = c(),
        upper = c(),
        optFn = function(par, NR){
            dnbinom(NR, mu = par[1], size = par[2])
        }
    ), 
    lossLateRep = list( # windows are sampled from two distributions at NR = NR_ploidy, NR_ploidy / 2
        par   = c(peakFraction = 0.8), # ~equivalent to replication modelType peakIsReplicated, or a number of windows with a single-copy loss
        lower = c(peakFraction = 0.05),
        upper = c(peakFraction = 0.95), # TODO: reduce upper to 0.9 if overfitting single peaks (also below)
        optFn = function(par, NR){
            dnbinom(NR, mu = par[1] / 2, size = par[2]) * (1 - par[3]) + 
            dnbinom(NR, mu = par[1],     size = par[2]) *      par[3]
        }
    ), 
    gain = list( # windows are sampled from two distributions at NR = NR_ploidy, NR_ploidy_CNC==+1
        par   = c(peakFraction = 0.8), # thus, cell has a number of windows with a single-copy gain
        lower = c(peakFraction = 0.05),
        upper = c(peakFraction = 0.95),
        optFn = function(par, NR){
            dnbinom(NR, mu = par[1],                                 size = par[2]) *      par[3] +             
            dnbinom(NR, mu = par[1] * (env$PLOIDY + 1) / env$PLOIDY, size = par[2]) * (1 - par[3]) 
        }
    ), 
    earlyRep = list( # windows are sampled from two distributions at NR = NR_ploidy, NR_ploidy * 2
        par   = c(peakFraction = 0.8), # ~equivalent to replication modelType peakIsUnreplicated
        lower = c(peakFraction = 0.05),
        upper = c(peakFraction = 0.95),
        optFn = function(par, NR){
            dnbinom(NR, mu = par[1],     size = par[2]) *      par[3] +             
            dnbinom(NR, mu = par[1] * 2, size = par[2]) * (1 - par[3])  
        }
    ),     
    lossAndGain = list( # windows are sampled from three distributions at NR = NR_ploidy_CNC==-1, NR_ploidy, NR_ploidy_CNC==+1
        par   = c(peakFraction = 0.6),   # thus, cell has a number of windows with each of a single-copy gain and single-copy loss
        lower = c(peakFraction = 0.05),  # for parsimonious fitting, force the loss and gains to the same overall frequency
        upper = c(peakFraction = 0.8),   # insist on at least 10% in CNV peaks to minimize symmetric overfitting of single peaks
        optFn = function(par, NR){
            dnbinom(NR, mu = par[1] * (env$PLOIDY - 1) / env$PLOIDY, size = par[2]) * ((1 - par[3]) / 2) +
            dnbinom(NR, mu = par[1],                                 size = par[2]) *       par[3] +             
            dnbinom(NR, mu = par[1] * (env$PLOIDY + 1) / env$PLOIDY, size = par[2]) * ((1 - par[3]) / 2)
        }
    )
) 
getCellDensity <- function(NR_wmsf){ # return the actual NR_wms frequency distribution, mixture models are fit to it
    density <- data.table(NR = NR_wmsf)[, .N, by = "NR"]
    NR_working <- with(density, { min(NR):max(NR) })
    density <- merge(data.table(NR = NR_working), density, all.x = TRUE)
    density[is.na(N), N := 0]
    density[, density := N / sum(N)]
    density
}
#=====================================================================================

#=====================================================================================
# functions that help adjust the window size to provide reliable copy number state shifts
#-------------------------------------------------------------------------------------
# determine the least number of reads per window for a cell's data to support a robust HMM
getMinWindowCount <- function(modal_NA){
    ploidyFactor <- (modal_NA + 0.5) / modal_NA
    (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# require that most windows have non-zero count, with allowance for rare homozygous losses and chrY
# similarly, the window median must have a usable count
# require a higher windowPower for a cell that cannot achieve these metrics
checkForExcessZeros <- function(cell_id, windowPower){
    windowSize <- 2 ** windowPower
    NR_w <- collapseVector(raw_counts[[cell_id]], windowSize)
    NR_wf <- NR_w[!is.na(NR_w)]
    mean(NR_wf < 1) > 0.2 || # return TRUE (reject) if too many zero-count windows
    median(NR_wf) < 5        # or if median is too low such that most windows have too few counts
}

# set the per-cell minimum window size as the number of bins needed to obtain a mean raw count >=minWindowCount
windowSizes <- 2 ** windowPowers
bins <- rowRanges$autosome & rowRanges$mappability >= env$MIN_MAPPABILITY
getMinWindowPower <- function(cell_id, minWindowCount){
    NR_b <- raw_counts[[cell_id]][bins]
    NR_b_mean <- mean(NR_b, na.rm = TRUE)
    if(is.na(NR_b_mean) || NR_b_mean == 0) return(env$MAX_WINDOW_POWER)
    windowSize <- ceiling(minWindowCount / NR_b_mean) # in number of bins, not bp
    if(is.na(windowSize) || windowSize > 2 ** env$MAX_WINDOW_POWER) return(env$MAX_WINDOW_POWER)
    windowSize <- windowSizes[min(which(windowSizes >= windowSize))] # thus, allow 2**(0:MAX_WINDOW_POWER) bins per window
    windowPower <- log2(windowSize)
    while(checkForExcessZeros(cell_id, windowPower) && 
          windowPower < env$MAX_WINDOW_POWER) windowPower <- windowPower + 1
    max(windowPower, env$MIN_WINDOW_POWER) # thus, allow 2**(MIN_WINDOW_POWER:MAX_WINDOW_POWER) bins per window
}

# use nrModelTypes, above, to determine if a windowPower results in an acceptable width of the NR_wm distribution
checkNA95 <- function(windowPower, NA95, scalar = NULL){
    if(is.null(scalar)) scalar <- if(windowPower == env$MAX_WINDOW_POWER) 1.5 else 1
    MAX_NA95 <- env$MAX_NA95 * scalar # be extra generous on the largest windowPower to allow more cells to be "not keep", rather than "bad"
    NA95 <= MAX_NA95
}
fitCellWindows <- function(windowPower, windows, NR_wm, shape = 1){
    NR_wmf <- excludeOutliers(NR_wm, min = 1) # thus, chrom-to-chrom, baseline-to-CNV, and rare bad windows won't have undue influence
    ER_modal_NA <- peakValue(NR_wmf) # ER_modal_NA may NOT be the same as ER_ploidy 
    sexChrom <- windows[, chrom %in% c("chrX", "chrY")] # use autosomes only to minimize the expected impact of sex chromosome "aneuploidy"
    RPA <- ER_modal_NA / env$PLOIDY 
    NR_min <- max(RPA / 10, 2) # suppress CN==0 bins from the fit
    NR_max <- RPA * 5          # suppress high CN bins from the fit, but allow 4N replication state
    NR_wms <- getNR_wms(NR_wm, shape)
    I <- which(!is.na(NR_wms) & NR_wms >= NR_min & NR_wms <= NR_max & !sexChrom)
    NR_wmsf <- round(NR_wms[I], 0) 
    density <- getCellDensity(NR_wmsf)
    fits <- do.call(rbind, lapply(names(nrModelTypes), function(modelTypeName){
        modelType <- nrModelTypes[[modelTypeName]]
        fit <- optim( # by fitting both ER_unreplicated and theta, we arrive at best mu and size parameters for nbinom for each model type
            c(
                ER_modal_NA = ER_modal_NA, 
                theta = 30, 
                modelType$par
            ),
            function(par){
                d <- modelType$optFn(par, density$NR) 
                sqrt(mean((d - density$density) ** 2))
            },
            method = "L-BFGS-B",
            lower = c(
                ER_modal_NA = ER_modal_NA * 0.9, # constrain the allowed movement of mu from the mode, i.e., peak
                theta = 1, 
                modelType$lower
            ),
            upper = c(
                ER_modal_NA = ER_modal_NA * 1.1, 
                theta = 1000, 
                modelType$upper
            )
        )    
        ER_modal_NA  <- fit$par[1]
        theta        <- fit$par[2]
        peakFraction <- if(modelTypeName == "singlePeak") 1 else fit$par[3]
        q95_modal_NA <- qnbinom(0.95, mu = ER_modal_NA, size = theta)
        RPA <- ER_modal_NA / env$PLOIDY
        NA95 <- (q95_modal_NA - ER_modal_NA) / RPA # how far the 95% CI is from the peak, in number of alleles
        data.table( # one of these is ultimately chosen to become the working cell object
            ER_modal_NA   = ER_modal_NA,
            theta         = theta,
            peakFraction  = peakFraction,            
            nrModelType   = modelTypeName,
            ER_ploidy     = ER_modal_NA,
            RPA           = RPA, # later we may learn RPA is 2-fold off for replicating cells
            rmsd          = fit$value,
            NA95          = NA95,
            allow         = checkNA95(windowPower, NA95)
        )
    }))
    bestFit <- as.list(fits[which.min(rmsd)[1]])
    bestFit$shape  <- shape  # on first pass, shape is 1 and NR_wms == NR_wm
    bestFit$NR_wms <- NR_wms # on second pass, chromosomes are shaped, NR_wms reflects that but NOT GC bias correction
    bestFit$density <- density
    bestFit
}
getCellWindows <- function(cell){ # parse a cell's windows at a given windowPower
    windowSize <- 2 ** cell$windowPower
    windows <- windows[[cell$windowPower + 1]]
    mappability <- windows[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_wm <- unname(unlist(sapply(constants$chrom, function(chrom){
        collapseVector(
            raw_counts[[cell$cell_id]][rowRanges$chrom == chrom], 
            windowSize
        ) / mappability[windows$chrom == chrom]
    })))
    fitCellWindows(cell$windowPower, windows, NR_wm)
}
plotCellWindows <- function(cell) {
    plotCountDensity(cell)
    plotWindows_counts(cell, paste0("windowPower", cell$windowPower), 
                       cell$windows$unshaped$NR_wms, increment = FALSE, col = "useGcColor") 
}
setCellWindows <- function(cell){ # establish the optimal window power for a cell
    workingStep <<- "setCellWindows"
    cell$windowPower <- cell$minWindowPower
    cell$windows <- list()
    cell$windows$unshaped <- getCellWindows(cell)
    plotCellWindows(cell)
    while(!cell$windows$unshaped$allow && # stop when we get a nice, appropriately tight distribution
          cell$windowPower < env$MAX_WINDOW_POWER){
        cell$windowPower <- cell$windowPower + 1
        cell$windows$unshaped <- getCellWindows(cell)
        plotCellWindows(cell)
    }
    if(cell$windowPower < env$MAX_WINDOW_POWER){ # also make input plots for all windowPowers up to one ABOVE the power used for fitting
        tmp <- copy(cell)
        tmp$windowPower <- tmp$windowPower + 1
        tmp$windows$unshaped <- getCellWindows(tmp)
        plotCellWindows(tmp)
        cell$nextWindows <- tmp$windows$unshaped # and save the lower resolution windows for possible use by late replicating cells
    }
    plotNumber <<- plotNumber + 1
    if(cell$windows$unshaped$allow) {
        setChromCN(cell) # pass the baton
    } else {
        cell$badCell <- TRUE # unusable data, abort cell with no further analysis
        cell$nextWindows <- NULL
        plotWindows_counts(cell, "input", cell$windows$unshaped$NR_wms, col = "badCell")
        cell
    }
}
#=====================================================================================

#=====================================================================================
# shared cell fitting functions, for GC bias effects and copy number by HMM
#-------------------------------------------------------------------------------------
getRepModelKey <- function(composite) if(composite) "composite" else "sequential"
getGcFit <- function(cell, composite) {
    if(composite) cell$replicationModel[[shapeKey]]$gc_fit
    else cell$windows[[shapeKey]]$sequential$gc_fit # this fit squashes replication GC effects
}
fitGcBias <- function(cell, name, NR_wm, windowI = TRUE, chroms = NULL, 
                      windowCN = NULL, col_w = NULL, composite = FALSE){
    windows <- windows[[cell$windowPower + 1]]
    if(!is.null(chroms)) windowI <- windowI & windows[, chrom %in% chroms]        
    gc_wf <- windows[windowI, gc_fraction]
    NR_wmf <- NR_wm[windowI]
    if(is.null(windowCN)) windowCN <- cell$ploidy
    if(length(windowCN) > 1 && length(windowCN) != length(NR_wmf)) windowCN <- windowCN[windowI]
    repKey <- getRepModelKey(composite)
    cell$windows[[shapeKey]][[repKey]]$gc_fit <- new_nbinomCountsGC2(NR_wmf, gc_wf, binCN = windowCN)
    plotWindows_counts(cell, name, NR_wm, col = col_w)   
    plotGcBias(cell, name, composite, gc_wf, NR_wmf, windowCN = windowCN)
    cell
}
solveSquashedHMM <- function(cell, name, NR_wm, windowI = TRUE){ # i.e., the CN HMM that is the first of the sequential model
    windows <- windows[[cell$windowPower + 1]]
    gc_wf <- windows[windowI, gc_fraction]
    NR_wmf <- NR_wm[windowI]
    cell$windows[[shapeKey]]$sequential$HMM <- viterbi(
        cell$windows[[shapeKey]]$sequential$gc_fit, 
        NR_wmf, # note, could simply be NR_wm if windowI == TRUE
        gc_wf, 
        maxCN = maxSequentialCN, 
        transProb = 1e-8,
        chroms = windows[windowI == TRUE, chrom], 
        asRle = FALSE
    )$cn
    plotWindows_cn(cell, name, FALSE, gc_wf, NR_wmf, windowI, cell$windows[[shapeKey]]$sequential$HMM)
    cell
}
#=====================================================================================

#=====================================================================================
# use a replication-unaware algorithm that squashes the replication GC bias effect
# to obtain an initial copy number estimates across the genome, i.e., CN, not NA
#-------------------------------------------------------------------------------------
setChromCN <- function(cell){ 
    workingStep <<- "setChromCN"
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # determine a first crude aneuploidy copy number estimate for each whole chromosome
    # this method based on the Mann-Whitney test does not assume anything about the NR_wm distribution
    # and is therefore resistant to different replication states, i.e., fractionS
    NR_wmsf <- excludeOutliers(cell$windows[[shapeKey]]$NR_wms, min = 2) # could simply be NR_wm if not shaped yet
    if(length(NR_wmsf) > env$MAX_FIT_WINDOWS) NR_wmsf <- sample(NR_wmsf, env$MAX_FIT_WINDOWS)
    RPA_wmsf <- NR_wmsf / cell$ploidy
    NR_wmsf_cn <- lapply(chromCNs, "*", RPA_wmsf)
    cell$tmp$CN_c <- lapply(chroms, function(chrom){
        NR_wmscl <- excludeOutliers(cell$windows[[shapeKey]]$NR_wms[windows$chrom == chrom], min = 2)
        if(length(NR_wmsf_cn) <= 1 || length(NR_wmscl) == 0 || median(NR_wmscl) <= 1) return(0)
        p <- sapply(NR_wmsf_cn, function(NR_g) wilcox.test(NR_g, NR_wmscl, exact = FALSE)$p.value)
        chromCNs[which.max(p)] # the copy number that gives the best distribution match between chrom and genome
    })
    names(cell$tmp$CN_c) <- chroms

    # use the crude estimates to select which chromosomes are used for performing a GC bias fit and CN HMM
    # only use chromosomes close to the aggregate distribution, i.e., those likely to be mostly CN==ploidy
    # unused chromosome windows are colored red in the QC plot
    chromFitI <- unlist(cell$tmp$CN_c) == cell$ploidy
    col_w <- sapply(windows$chrom, function(chrom) {
        if(cell$tmp$CN_c[[chrom]] == cell$ploidy) defaultPointColor else rgb(1, 0, 0, 0.1)
    })    # vvvvvvv this GC bias fit squashes the replication effect on purpose
    cell <- fitGcBias(cell, "input", cell$windows[[shapeKey]]$NR_wms, chroms = chroms[chromFitI], col_w = col_w, composite = FALSE)
    cell <- solveSquashedHMM(cell, "squashed", cell$windows[[shapeKey]]$NR_wms)

    # revise the gc_fit and HMM in one iteration cycle
    name <- "revised"
    cell <- fitGcBias(cell, name, cell$windows[[shapeKey]]$NR_wms, 
                       windowCN = cell$windows[[shapeKey]]$sequential$HMM, col_w = col_w, composite = FALSE)
    cell <- solveSquashedHMM(cell, name, cell$windows[[shapeKey]]$NR_wms)
    setReplicationModel(cell) # pass the baton
}
#=====================================================================================

#=====================================================================================
# establish the common parameters of a set of replication models to be applied to individual cells
# these models are used on NR_wm rescaled by inferred CN and thus isolate the replication effect
#-------------------------------------------------------------------------------------
modelTypes <- list(
    notReplicating     = 0, # allowable fractionsS values by peak type
    peakIsUnreplicated = seq(0.05, 0.75, 0.025),
    peakIsReplicated   = seq(0.25, 0.95, 0.025)
)
modelValues <- do.call(rbind, lapply(names(modelTypes), function(modelType){
    data.table( # pre-assemble a table of all models
        modelType = modelType,
        fractionS = modelTypes[[modelType]]
    )
}))
modelTypeColors <- list( # for plotting
    notReplicating      = "blue",
    peakIsUnreplicated  = "green",
    peakIsReplicated    = "red"
)
fractionSFile <- file.path(env$ACTION_DIR, "fractionS_table.rds")
fractionSLookup <- readRDS(fractionSFile) # for converting overall fractionS to fraction of 0, 1, and 2 alleles replicated
#=====================================================================================

#=====================================================================================
# support functions for fitting replication models
#-------------------------------------------------------------------------------------
sumLogLikelihoods <- function(x){
    x[is.na(x)] <- 0
    x[x == -Inf] <- -800 # suppress outliers      
    sum(x)
}
getModelER_unreplicated <- function(modelType, cell){ # set read count expectations (ER) for a combination of modelType and cell$ER_modal_NA
    with(cell$windows[[shapeKey]], { switch(
        modelType,
        notReplicating     = ER_modal_NA,
        peakIsUnreplicated = ER_modal_NA,    # the modal window count corresponds to unreplicated DNA
        peakIsReplicated   = ER_modal_NA / 2 # the modal window count corresponds to   replicated DNA
    )})
}
fitReplicationModel <- function(modelType, fractionS, cell){
    # use a rich subset of the data to fit the candidate model to the data as best as possible
    # models are defined by replication peak type and fractionS increments    
    ER_unreplicated <- getModelER_unreplicated(modelType, cell)
    theta <- cell$windows[[shapeKey]]$theta
    weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, fractionS)
    optFn <- function(par) with(cell$tmp$density, { 
        d <- dnbinom(NR, mu = par[1] * 1,   size = par[2]) * weights[1] +   
             dnbinom(NR, mu = par[1] * 1.5, size = par[2]) * weights[2] + 
             dnbinom(NR, mu = par[1] * 2,   size = par[2]) * weights[3]  
        sqrt(mean((d - density) ** 2))
    })
    opt <- optim( # by fitting both ER_unreplicated and theta, we arrive at mu and size parameters for nbinom
        c(ER_unreplicated = ER_unreplicated, theta = theta), # adjusted from our initial estimate made from the peak, i.e., mode
        optFn,
        method = "L-BFGS-B",
        lower = c(ER_unreplicated = ER_unreplicated * 0.9, theta = 1),
        upper = c(ER_unreplicated = ER_unreplicated * 1.1, theta = theta * 10)
    )
    opt$par[3] <- opt$value
    opt$par
}
fillAllReplicationModels <- function(cell){  # fit all required models
    models <- copy(modelValues)
    # if(!is.null(cell$replicationModel)){ # when refitting for ER_modal_NA and theta, just use the single established model
    #     models <- with(cell, { models[ modelType == replicationModel$modelType & fractionS == replicationModel$fractionS ] })
    # }
    models <- cbind(models, t(models[, mapply(function(modelType, fractionS){
        fitReplicationModel(modelType, fractionS, cell)
    }, modelType, fractionS)]))
    setnames(models, c("modelType", "fractionS", "ER_unreplicated", "theta", "rmsd"))
    models
}
getRepEmissProbs <- function(cell, model){ # used to solve an inital replication HMM
    ER_unreplicated <- unlist(model$ER_unreplicated)
    with(cell$tmp, { log(as.matrix(data.table( # here, model$ER_unreplicated has been fitted to reflect mu, not peak/mode
        dnbinom(NR_wmshf, mu = ER_unreplicated,       size = model$theta), 
        dnbinom(NR_wmshf, mu = ER_unreplicated * 1.5, size = model$theta), 
        dnbinom(NR_wmshf, mu = ER_unreplicated * 2,   size = model$theta)
    ))) })
}
setBestReplicationModel <- function(cell, models, name, windows, fitI){  # select the best model from among those still allowed
    bestModelByType <- models[, { 
        i <- which.min(rmsd)[1]
        .SD[i]
    }, by = "modelType"]

    # update the likelihood of each candidate model by fitting replication HMMs
    gc_wf <- windows[fitI, gc_fraction]
    bestModelByType <- merge(bestModelByType, bestModelByType[, {
        model <- .SD
        if(modelType == "notReplicating") { # always allow the default nonReplicating model
            gc_fit <- getGcFit(cell, composite = FALSE)
            RPA_wf <- predict(gc_fit, gc_wf)
            ER_whf <- RPA_wf * cell$ploidy
            gcModel <- optimize(function(theta){
                sumLogLikelihoods(log(dnbinom(cell$tmp$NR_wmshf, mu = ER_whf, size = theta)))
            }, c(1, model$theta * 10), maximum = TRUE)
            list(
                NAR = list(rep(0, length(fitI))),
                logLikelihood = gcModel$objective,
                theta2 = gcModel$maximum,
                gc_fit = list(gc_fit),
                gcRatio = as.double(NA),
                fractionS2 = 0.0,
                fractionFullyReplicated = 1.0, 
                peakRatio = 0.0,
                allow = TRUE
            )         
        } else { # always solve the HMM for the best candidate early and late replicating cell models

            # run a first replication HMM to get a first estimate of the window replication state
            emissProbs <- getRepEmissProbs(cell, model)
            hmm <- new_hmmEPTable(emissProbs, transProb = repTransProb, keys = windows$chrom[fitI])
            nar_hmm <- keyedViterbi(hmm) - 1 # based on NR_wmshf
            nar_hmm[is.na(nar_hmm)] <- 0 # happens mainly on chrY when there is only 1 window to analyze

            # update replication HMM with both replication and systematic GC bias components
            na_hmm <- cell$ploidy + nar_hmm            
            gc_fit <- new_nbinomCountsGC2(cell$tmp$NR_wmshf, gc_wf, binCN = na_hmm)
            RPA_wf <- predict(gc_fit, gc_wf)
            model$ER_unreplicated <- list(RPA_wf * cell$ploidy)
            ER_whf <- RPA_wf * na_hmm
            model$theta <- optimize(function(theta){
                sumLogLikelihoods(log(dnbinom(cell$tmp$NR_wmshf, mu = ER_whf, size = theta)))
            }, c(1, model$theta * 10), maximum = TRUE)$maximum
            emissProbs <- getRepEmissProbs(cell, model)
            hmm <- new_hmmEPTable(emissProbs, transProb = repTransProb, keys = windows$chrom[fitI])
            nar_hmm <- keyedViterbi(hmm) - 1 # based on NR_wmshf
            nar_hmm[is.na(nar_hmm)] <- 0

            # require a user-specified minimal fractionS
            gcRatio    <- as.double(NA)
            fractionS2 <- as.double(NA)
            peakRatio  <- as.double(NA)            
            nar_rescaled <- nar_hmm * cell$windows[[shapeKey]]$sequential$HMM[fitI] / cell$ploidy
            fractionS2 <- sum(nar_rescaled, na.rm = TRUE) / (sum(cell$windows[[shapeKey]]$sequential$HMM[fitI], na.rm = TRUE))
            minorFraction <- 0.5 - abs(fractionS2 - 0.5)   

            # further require substantial 4N content in early replicating cells (avoids overcalling of 3N as NAR==1)
            fractionFullyReplicated <- if(modelType == "peakIsUnreplicated"){
                sum(nar_hmm == 2, na.rm = TRUE) / sum(fitI, na.rm = TRUE)              
            } else 1 - minorFraction

            # require that replicated DNA have a higher GC content than unreplicated DNA
            windowIsReplicated <- nar_hmm > 0              
            gc_replicated   <- windows[fitI][windowIsReplicated == TRUE,  mean(gc_fraction, na.rm = TRUE)]
            gc_unreplicated <- windows[fitI][windowIsReplicated == FALSE, mean(gc_fraction, na.rm = TRUE)]
            gcRatio <- gc_replicated / gc_unreplicated

            # require that the windows called as fully replicated have the expect copy number
            NR_nar_4N <- cell$tmp$NR_wmshf[nar_hmm == 2] # fully replicated
            NR_nar_2N <- cell$tmp$NR_wmshf[nar_hmm == 0] # unreplicated
            peakValue2N <- if(length(NR_nar_2N) < 100) median(NR_nar_2N, na.rm = TRUE) else peakValue(NR_nar_2N)
            peakValue4N <- if(length(NR_nar_4N) < 100) median(NR_nar_4N, na.rm = TRUE) else peakValue(NR_nar_4N)
            peakRatio   <- peakValue4N / peakValue2N # expected value == 2

            # check if the replication model meets expectations
            if(minorFraction < env$MIN_FRACTION_S || 
               fractionFullyReplicated < 0.01 ||
               is.na(gcRatio) || gcRatio < 1.05 || 
               is.na(peakRatio) || peakRatio < 1.75){
                list( 
                    NAR = list(),
                    logLikelihood = as.double(NA),
                    theta2 = as.double(NA),
                    gc_fit = list(),
                    gcRatio = gcRatio,
                    fractionS2 = fractionS2,
                    fractionFullyReplicated = fractionFullyReplicated,
                    peakRatio = peakRatio,
                    allow = FALSE
                )

            # get the final best fit of the replication model (it may still be rejected if nonReplicating is a better fit)    
            } else {
                na_hmm <- cell$ploidy + nar_hmm            
                gc_fit <- new_nbinomCountsGC2(cell$tmp$NR_wmshf, gc_wf, binCN = na_hmm)
                RPA_wf <- predict(gc_fit, gc_wf)
                ER_whf <- RPA_wf * na_hmm
                gcModel <- optimize(function(theta){
                    sumLogLikelihoods(log(dnbinom(cell$tmp$NR_wmshf, mu = ER_whf, size = theta)))
                }, c(1, model$theta * 10), maximum = TRUE)
                list(
                    NAR = list({
                        NAR <- rep(0, length(fitI))
                        NAR[fitI] <- nar_hmm
                        round(NAR * cell$windows[[shapeKey]]$sequential$HMM / cell$ploidy, 0) # rescale NAR to hmm
                    }),
                    logLikelihood = gcModel$objective,
                    theta2 = gcModel$maximum,
                    gc_fit = list(gc_fit),
                    gcRatio = gcRatio,
                    fractionS2 = fractionS2,
                    fractionFullyReplicated = fractionFullyReplicated,
                    peakRatio = peakRatio,
                    allow = TRUE
                )
            }
        }        
    }, by = "modelType"])

    # select and characterize the best model
    cell$replicationModel[[shapeKey]] <- as.list(bestModelByType[allow == TRUE][which.max(logLikelihood)[1]])
    cell$windows[[shapeKey]]$sequential$NAR <- unlist(cell$replicationModel[[shapeKey]]$NAR)
    cell$windows[[shapeKey]]$sequential$fractionS <- cell$replicationModel[[shapeKey]]$fractionS2 # as obtained from the sequential HMM >> NAR fits
    cell$replicationModel[[shapeKey]] <- with(cell$replicationModel[[shapeKey]], list(
        modelType = modelType,
        fractionS = fractionS, # as obtained from the initial mixture model with stepped fractionS values
        ER_unreplicated = ER_unreplicated, 
        theta = theta2,
        gc_fit = gc_fit[[1]],
        logLikelihood = logLikelihood,
        gcRatio = gcRatio,
        peakRatio = peakRatio,
        cellIsReplicating = modelType != "notReplicating"
    ))
    if(cell$replicationModel[[shapeKey]]$modelType == "peakIsReplicated"){ # adjust values now that we now that peak == CN4, not CN2
        cell$modal_NA <- cell$ploidy * 2
        cell$windows[[shapeKey]]$ER_ploidy <- cell$windows[[shapeKey]]$ER_modal_NA / 2
        cell$windows[[shapeKey]]$RPA <- cell$windows[[shapeKey]]$ER_modal_NA / cell$modal_NA 
        cell$windows[[shapeKey]]$NA95 <- cell$windows[[shapeKey]]$NA95 * 2
        cell$windows[[shapeKey]]$allow <- with(cell, { checkNA95(windowPower, windows[[shapeKey]]$NA95, scalar = 1.25) })        
    }
    plotCellFit(cell, name, cell$replicationModel[[shapeKey]], bestModelByType) # based on NR_wmshfl
    plotNumber <<- plotNumber + 1
    plotGcBias(cell, name, FALSE, gc_wf, cell$windows[[shapeKey]]$NR_wms[fitI], 
               windowCN = with(cell$windows[[shapeKey]]$sequential, {HMM[fitI] + NAR[fitI]}), 
               gc_fit = cell$replicationModel[[shapeKey]]$gc_fit, windowI = fitI) 
    plotWindows_cn(cell, name, FALSE, gc_wf, cell$windows[[shapeKey]]$NR_wms[fitI], windowI = fitI, 
                   hmm = cell$windows[[shapeKey]]$sequential$HMM[fitI], gc_fit = cell$replicationModel[[shapeKey]]$gc_fit)
    cell
}
#=====================================================================================

#=====================================================================================
# convert a cell's HMM + NAR results to derivative values
#-------------------------------------------------------------------------------------
setHMMProfiles <- function(cell, composite){
    repKey <- getRepModelKey(composite)
    informative <- with(cell$windows[[shapeKey]][[repKey]], { HMM > 0 & HMM < maxCompositeCN }) # the HMM CN states that can be used and trusted in aggregated calculations  
    cell$windows[[shapeKey]][[repKey]]$fractionS <- with(cell$windows[[shapeKey]][[repKey]], { 
        sum(NAR[informative], na.rm = TRUE) / 
        sum(HMM[informative], na.rm =TRUE) 
    }) # FAR for the whole cell 
    NA_ <- with(cell$windows[[shapeKey]][[repKey]], { HMM + NAR })
    NR_wmsh <- with(cell, { windows[[shapeKey]]$NR_wms * ploidy / windows[[shapeKey]][[repKey]]$HMM }) # correct NR toward ploidy, for analyzing replication
    NR_wmsr <- with(cell, { windows[[shapeKey]]$NR_wms * windows[[shapeKey]][[repKey]]$HMM / NA_ })    # correct NR toward unreplicated, for analyzing CNVs
    ploidy <- cell$ploidy  
    windows <- windows[[cell$windowPower + 1]]
    gc_w <- windows[, gc_fraction]
    gc_fit <- cell$replicationModel[[shapeKey]]$gc_fit # the truest available RPA estimate
    RPA_peak <- predict(gc_fit, gc_w, type = 'adjustedPeak')
    cell$windows[[shapeKey]][[repKey]]$CN  <-  NR_wmsr / RPA_peak # calculated copy number after accounting for replication
    cell$windows[[shapeKey]][[repKey]]$FAR <- (NR_wmsh / RPA_peak - ploidy) / ploidy # fraction replicated per window
    plotWindows_CN_isolated( cell, paste(repKey, "isolated"), composite)
    plotWindows_FAR_isolated(cell, paste(repKey, "isolated"), composite)

    # TODO: commit/save HMM + NAR profiles of replicating cells, or, even better, its allele votes 

    # calculate a post-hoc assessment of the quality of the fully fitted copy number model
    x <- with(cell, { data.table(NR = windows[[shapeKey]]$NR_wms, RPA = RPA_peak, NA_ = NA_, HMM = windows[[shapeKey]][[repKey]]$HMM) })
    x <- x[HMM %in% 1:3, .(
        metric = if(.N >= 2) sd(NR / RPA / NA_, na.rm = TRUE) * sqrt(NA_[1]) else 0.5,
        weight = .N
    ), by = NA_]
    cell$windows[[shapeKey]][[repKey]]$cnsd <- weighted.mean(x$metric, x$weight, na.rm = TRUE)
    cell$windows[[shapeKey]][[repKey]]$keep <- cell$windows[[shapeKey]][[repKey]]$cnsd <= env$KEEP_THRESHOLD # thus, two SD = 0.8 CN unit when KEEP_THRESHOLD == 0.4
    cell
}
#=====================================================================================

#=====================================================================================
# find and characterize a cell's HMM segments
#-------------------------------------------------------------------------------------
replicatingStates <- fread(file.path(env$ACTION_DIR, "composite_states.csv"))
nonReplicatingStates <- replicatingStates[NAR == 0]
getHmmSegments <- function(cell, composite){
    repKey <- getRepModelKey(composite)
    cellStates <- copy(cell$states[[shapeKey]])
    setkey(cellStates, "NAME")
    windowSizeBases <- 2 ** cell$windowPower * 2e4 # bin size fixed at 20kb
    if(nrow(windows[[cell$windowPower + 1]]) != length(cell$windows[[shapeKey]][[repKey]]$stateI))
        stop("getHmmSegments error: windows and stateI lengths don't match")
    if(any(is.na(cell$windows[[shapeKey]][[repKey]]$stateI)))
        stop("getHmmSegments error: missing values for stateI")
    segments <- windows[[cell$windowPower + 1]][, { # segments are called per-chromosome, i.e, they obey chromosome ends
        segments <- with(cell, { rle(cellStates[windows[[shapeKey]][[repKey]]$stateI[.I], NAME]) })
        x <- cellStates[segments$values]
        ends <- cumsum(segments$lengths)
        starts <- head(c(1, ends + 1), -1)
        x[, ":="(
            chromStartI = starts,
            chromEndI = ends,
            start = start[starts],
            end = end[ends]
        )]
        x
    }, by = chrom]
    segments[, ":="(
        FAR = ifelse(CN == 0, 0, NAR / CN),
        CNC_P1 = NULL, # clean up unneeded columns
        CNC_P2 = NULL
    )]
    x <- list( # return a list with three different kinds of segments
        N_ALLELE = segments,     # i.e., all HMM segments regardless of the reason for each N_ALLELE break
        CN = segments[, .SD[, .( # i.e., all HMM breaks that resulted from a CN switch
            NAME = list(sort(unique(NAME))),
            CN   = CN[1],
            NAR  = list(sort(unique(NAR))),
            FAR  = list(sort(unique(FAR))),
            N_ALLELES = list(sort(unique(N_ALLELES))),
            chromStartI = min(chromStartI),
            chromEndI = max(chromEndI),
            start = min(start),
            end = max(end)
        ), by = rleidv(CN)], by = "chrom"],
        FAR = segments[, .SD[, .( # i.e., all HMM breaks that resulted from an altered replication state
            NAME = list(sort(unique(NAME))),
            CN   = list(sort(unique(CN))),
            NAR  = list(sort(unique(NAR))),
            FAR  = FAR[1],
            N_ALLELES = list(sort(unique(N_ALLELES))),
            chromStartI = min(chromStartI),
            chromEndI = max(chromEndI),
            start = min(start),
            end = max(end)
        ), by = rleidv(FAR)], by = "chrom"] # rleidv() creates an ID column based on a vector rle call
    )
    x$N_ALLELE[, ID := .I]
    setnames(x$CN, 2, "ID")
    setnames(x$FAR, 2, "ID")
    x
}
#=====================================================================================

#=====================================================================================
# fit a cell's data to establish the initial parameters of it's replication model
#-------------------------------------------------------------------------------------
setReplicationModel <- function(cell){
    workingStep <<- "setReplicationModel"
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # use CNV-corrected counts to establish cell-specific likelihoods for replication models
    cell$tmp$NR_wmsh <- with(cell, { windows[[shapeKey]]$NR_wms * ploidy / windows[[shapeKey]]$sequential$HMM }) 
    fitI <- with(cell, { windows[[shapeKey]]$sequential$HMM > 0 & 
                         windows[[shapeKey]]$sequential$HMM < maxSequentialCN & 
                         !is.na(tmp$NR_wmsh) }) 
    cell$tmp$NR_wmshf <- round(cell$tmp$NR_wmsh[fitI], 0)
    ER_modal_NA <- cell$windows[[shapeKey]]$ER_modal_NA
    cell$tmp$NR_wmshfl <- with(cell$tmp, { NR_wmshf[NR_wmshf >= 2 & # suppress CN == 0
                                                    NR_wmshf >= ER_modal_NA / 2 / 4 & # limit outlier effects
                                                    NR_wmshf <= ER_modal_NA * 2 * 2] })
    cell$tmp$density <- getCellDensity(cell$tmp$NR_wmshfl)

    # optimize theta for three distinct types of models, at multiple fixed fractionS values
    # and ER_(un)replicated values determine by model type
    # pick the best model type, with its best parameters
    isFirstPass <- length(cell$windows[[shapeKey]]$shape) == 1
    name <- if(isFirstPass) "initial" else "shaped"    
    models <- fillAllReplicationModels(cell)
    cell <- setBestReplicationModel(cell, models, name, windows, fitI)
    cellIsReplicating <- cell$replicationModel[[shapeKey]]$cellIsReplicating

    # re-check fit resolution for late replicating cells where we adjusted the peak from CN2 to CN4
    # if too wide, and if possible, move to the next highest window power that we cached above
    # if(cell$cellIsReplicating && !cell$windows[[shapeKey]]$allow && !is.null(cell$nextWindows)){
    #     cell$windowPower <- cell$windowPower + 1
    #     cell$windows[[shapeKey]] <- cell$nextWindows
    #     cell$nextWindows <- NULL # thus, we will only do this once
    #     cell$replicationModel <- NULL
    #     cell$tmp <- NULL
    #     return(setChromCN(cell))
    # }

    # finalize the sequential model and solve the composite HMM for replicating cells
    cell$windows$fitI <- fitI    
    cell <- setHMMProfiles(cell, composite = FALSE)
    cell$states[[shapeKey]] <- copy(if(cellIsReplicating) replicatingStates else nonReplicatingStates)
    setkey(cell$states[[shapeKey]], "NAME")
    cell$states[[shapeKey]][, i := .I] # this construction mimics what the composite HMM will do, below
    stateNames <- with(cell$windows[[shapeKey]]$sequential, { paste0("CN", HMM, "_NAR", NAR) })
    cell$windows[[shapeKey]]$sequential$stateI <- cell$states[[shapeKey]][stateNames, i]
    cell$segments[[shapeKey]]$sequential <- getHmmSegments(cell, FALSE)
    if(cellIsReplicating) solveCompositeHmm(cell) else finishShapeModel(cell)
}
#=====================================================================================

#=====================================================================================
# support functions for replication state profiles as a function of fractionS + window GC
# i.e., that add P(GC | NAR) emission probability component to the composite HMM for replicating cells
# derived from the intial fit of each individual cell
#-------------------------------------------------------------------------------------
nGcSteps <- 50
gcIndices  <- 1:nGcSteps
minPNar <- 0.01
getP_gc_nar <- function(cell, gc_w, iteration) with(cell, {
    repKey <- if(iteration == 1) "sequential" else "composite"
    cnNeutralI <- windows[[shapeKey]][[repKey]]$HMM[windows$fitI] == ploidy
    dt <- dcast(
        data.table(
            GCI = round(gc_w[windows$fitI][cnNeutralI] * nGcSteps, 0), # vector of filtered window GC indices
            NAR = windows[[shapeKey]][[repKey]]$NAR[windows$fitI][cnNeutralI]
        ), 
        GCI ~ NAR, 
        fun.aggregate = length, 
        fill = 0, 
        drop = FALSE
    )
    dt$all <- rowSums(dt[, 2:ncol(dt)])
    if(ploidy == 2) dt$replicated <- dt$all - dt[["0"]] 
    for(j in 2:ncol(dt)) dt[[j]] <- pmax(minPNar, dt[[j]] / sum(dt[[j]]))
    get_by_col <- function(nar) {
        x <- dt[[as.character(nar)]]
        if(is.null(x)) minPNar else x
    }
    get_by_mix <- function(nar1, nar2, w1, w2){
        x <- nar1 * w1 + nar2 * w2 # weighted average of two distributions
        x / sum(x) # such that all emission probability distributions sum to 1
    }
    if(ploidy == 1){
        dt$CN1_NAR0 <- get_by_col(0) # CN1, same as the NAR model
        dt$CN1_NAR1 <- get_by_col(1) 
        #-------------
        dt$CN0_NAR0 <- dt$all # CN0, one state, replication not meaningful
        #-------------
        dt$CN2_NAR0 <- dt$CN1_NAR0 # other CN, interpolate intermediate values
        dt$CN2_NAR1 <- get_by_mix(dt$CN1_NAR0, dt$CN1_NAR1, 1/2, 1/2)
        dt$CN2_NAR2 <- dt$CN1_NAR1
        #-------------
        dt$CN3_NAR0 <- dt$CN1_NAR0
        dt$CN3_NAR1 <- get_by_mix(dt$CN1_NAR0, dt$CN1_NAR1, 2/3, 1/3)
        dt$CN3_NAR2 <- get_by_mix(dt$CN1_NAR0, dt$CN1_NAR1, 1/3, 2/3)
        dt$CN3_NAR3 <- dt$CN1_NAR1
        #-------------
        dt$CN4_NAR0 <- dt$CN1_NAR0
        dt$CN4_NAR1 <- get_by_mix(dt$CN1_NAR0, dt$CN1_NAR1, 3/4, 1/4)
        dt$CN4_NAR2 <- get_by_mix(dt$CN1_NAR0, dt$CN1_NAR1, 1/2, 1/2)
        dt$CN4_NAR3 <- get_by_mix(dt$CN1_NAR0, dt$CN1_NAR1, 1/4, 3/4)
        dt$CN4_NAR4 <- dt$CN1_NAR1
    } else { 
        dt$CN2_NAR0 <- get_by_col(0) # CN2, same as the NAR model
        dt$CN2_NAR2 <- get_by_col(2)   
        if(is.null(dt[["1"]])) dt[["1"]] <- get_by_mix(dt$CN2_NAR0, dt$CN2_NAR2, 1/2, 1/2)
        dt$CN2_NAR1 <- get_by_col(1)
        #-------------
        dt$CN0_NAR0 <- dt$all # CN0, one state, replication not meaningful
        #-------------
        dt$CN1_NAR0 <- dt$CN2_NAR0 # other CN, interpolate intermediate values
        dt$CN1_NAR1 <- dt$replicated
        #-------------
        dt$CN3_NAR0 <- dt$CN2_NAR0
        dt$CN3_NAR1 <- get_by_mix(dt$CN2_NAR0, dt$CN2_NAR1, 1/2, 1/2)
        dt$CN3_NAR2 <- get_by_mix(dt$CN2_NAR1, dt$CN2_NAR2, 1/4, 3/4)
        dt$CN3_NAR3 <- dt$CN2_NAR2
        #-------------
        dt$CN4_NAR0 <- dt$CN2_NAR0
        dt$CN4_NAR1 <- get_by_mix(dt$CN2_NAR0, dt$CN2_NAR1, 1/2, 1/2)
        dt$CN4_NAR2 <- dt$CN2_NAR1
        dt$CN4_NAR3 <- get_by_mix(dt$CN2_NAR1, dt$CN2_NAR2, 1/2, 1/2)
        dt$CN4_NAR4 <- dt$CN2_NAR2  
    }
    dt$key <- as.character(dt$GCI)
    setkey(dt, "key")
    plotP_gc_nar(cell, dt, iteration)
    dt
})
#=====================================================================================

#=====================================================================================
# use the replication model and reshaped GC bias to solve a multi-state HMM with CN and replication components
#   use cell$NR_wms to normalize to chromosome shape, i.e., amplification skew
#   use gc_fit to determine ER per state as gc_fit(GC) * (CN + NAR = NA)
#   use P_unreplicated_wf to weight windows for timed replication potential based on GC content
# none of these are deterministic for a window, they are only trends determined from aggregated data
# thus, the model output is not dictated or forced at any window by upstream preparative work
#-------------------------------------------------------------------------------------
crossCheckCnvs <- function(query, reference, ploidy){ # each a data.table of segments, nrow(query) == 1, nrow(reference) = any 
    reference <- reference[CN != ploidy]
    if(nrow(reference) == 0) return(0)
    overlaps <- reference[
        chrom == query$chrom & 
        chromStartI <= query$chromEndI & 
        query$chromStartI <= chromEndI & 
        sign(CN - ploidy) == sign(query$CN - ploidy)
    ]
    if(nrow(overlaps) == 0) return(0)    
    queryI <- query$chromStartI:query$chromEndI
    queryLength <- query$chromEndI - query$chromStartI + 1 # return the fraction of query also called as a CNV in reference
    sum(overlaps[, sum(chromStartI:chromEndI %in% queryI), by = ID][[2]]) / queryLength
}
suppressFalseSequentialCnvs <- function(cell) with(cell$segments[[shapeKey]], {
    segmentsWereMasked <- FALSE
    checkReplicatingSegment <- function(i) with(cell$segments[[shapeKey]], {
        foundInComposite <- crossCheckCnvs(sequential$CN[i], composite$CN, cell$ploidy) > 0
        if(foundInComposite) return(NULL)
        segChrom <- sequential$CN[i, chrom]
        chromIRange <- with(sequential, { c(CN[i, chromStartI], CN[i, chromEndI]) })
        isCNZero <- sequential$CN[i, CN == 0] # avoid subsequent divide by zero during sequential NAR re-fit
        overrideI <- windows[[cell$windowPower + 1]][, {
            if(chrom == segChrom) {
                x <- 1:.N %between% chromIRange 
                n <- sum(x) # always pass whole-chromosome aneuploidy or especially large CNVs
                if(isCNZero || n == .N || n > 250) rep(FALSE, .N) else x
            } else rep(FALSE, .N)
        }, by = chrom][[2]]
        if(sum(overrideI) > 0){
            cell$windows[[shapeKey]]$sequential$HMM[overrideI] <<- cell$ploidy
            segmentsWereMasked <<- TRUE
        }
    })    
    with(cell$segments[[shapeKey]]$sequential, { 
        for(i in 1:nrow(CN)){
            if(CN[i, CN != cell$ploidy]) checkReplicatingSegment(i)
        }
    })    
    if(!segmentsWereMasked) return(cell)
    cell$windows[[shapeKey]]$sequential$NAR <- NULL
    cell
})
suppressFalseCompositeCnvs <- function(cell){    
    cell$segments[[shapeKey]]$composite <- getHmmSegments(cell, TRUE)
    checkReplicatingSegment <- function(i, requiredStates, fallbackState) with(cell$segments[[shapeKey]], {
        foundInSequential <- crossCheckCnvs(composite$CN[i], sequential$CN, cell$ploidy) > 0
        hasRequiredStates <- composite$CN[i, any(requiredStates %in% unlist(NAME))]
        if(foundInSequential && hasRequiredStates) return(NULL)
        segChrom <- composite$CN[i, chrom]
        chromIRange <- with(composite, { c(CN[i, chromStartI], CN[i, chromEndI]) })
        overrideI <- windows[[cell$windowPower + 1]][, {
            if(chrom == segChrom) {
                x <- 1:.N %between% chromIRange 
                n <- sum(x) # always pass whole-chromosome aneuploidy or especially large CNVs, even if fully replicated
                if(foundInSequential && (n == .N || n > 250)) rep(FALSE, .N) else x
            } else rep(FALSE, .N)
        }, by = chrom][[2]]
        if(sum(overrideI) > 0){
            fallBackStateI <- which(cell$states[[shapeKey]]$NAME == fallbackState)
            cell$windows[[shapeKey]]$composite$stateI[overrideI] <<- fallBackStateI
        }
    })    
    with(cell$segments[[shapeKey]]$composite, { if(cell$ploidy == 1){
        for(i in 1:nrow(CN))
            if(CN[i, CN == 2]) 
                checkReplicatingSegment(i, c("CN2_NAR1", "CN2_NAR2"), "CN1_NAR1") # i.e., require states not consistent with CN1
    } else if(cell$ploidy == 2){
        for(i in 1:nrow(CN))
            if(CN[i, CN == 1]) 
                checkReplicatingSegment(i, c("CN1_NAR0"), "CN2_NAR0") # i.e., require states not consistent with CN2
            else if(CN[i, CN == 3]) 
                checkReplicatingSegment(i, c("CN3_NAR2", "CN3_NAR3"), "CN2_NAR1")
            else if(CN[i, CN == 4]) 
                checkReplicatingSegment(i, c("CN4_NAR1", "CN4_NAR2", "CN4_NAR3", "CN4_NAR4"), "CN2_NAR2") 
    } })
    cell$segments[[shapeKey]]$composite <- getHmmSegments(cell, TRUE) # redo segments calls in case any CNV states were overridden
    cell
}
solveCompositeHmm <- function(cell){
    workingStep <<- "solveCompositeHmm"

    # get window values
    windows <- windows[[cell$windowPower + 1]]
    gc_w <- windows[, gc_fraction]
    gci_w <- as.character(round(gc_w * nGcSteps, 0))

    # solve two HMM cycles, using the early HMMs as guidance for replicating cells
    nStates <- nrow(cell$states[[shapeKey]]) - 1 # since we will remove CN5_NAR0 prior to composite HMM    
    theta <- cell$replicationModel[[shapeKey]]$theta
    RPA_mu <- predict(cell$replicationModel[[shapeKey]]$gc_fit, gc_w)
    NR_wms <- round(cell$windows[[shapeKey]]$NR_wms, 0)
    iterateCompositeHmm <- function(iteration, transProb){

        # adjust the states and state transition weights for this HMM iteration
        transProbs <- matrix(NA, nrow = nStates, ncol = nStates)
        for(fromI in 1:nStates) for(toI in 1:nStates){
            fromState <- cell$states[[shapeKey]][fromI]
            toState   <- cell$states[[shapeKey]][toI]
            transProbs[fromI, toI] <- if(fromI == toI) 1 - nStates * transProb
            else if(fromState$CN != toState$CN) {
                if(fromState$N_ALLELES != toState$N_ALLELES) transProb  # discourage CNVs
                else 1e-12 # STRONGLY discourage changing CN and NAR simultaneously to yield the same N_ALLELES
            } else repTransProb # changing NAR is much more likely in replicating cells than introducing a CNV
        }

        # establish the cell specific values for P(GC | NAR)
        P_gc_nar <- getP_gc_nar(cell, gc_w, iteration)
        P_gc_nar_w <- P_gc_nar[gci_w]   

        # calculate the emission probabilities using all prior fitting as guidance
        emissProbs <- sapply(cell$states[[shapeKey]][NAME != "CN5_NAR0"]$NAME, function(stateName){ # NAME	CN	NAR	N_ALLELES	CNC_P1	CNC_P2
            state <- cell$states[[shapeKey]][NAME == stateName]
            ER <- RPA_mu * if(state$N_ALLELES == 0) 0.05 else state$N_ALLELES
            pNA  <- dnbinom(NR_wms, mu = ER, size = theta)
            pNAR <- P_gc_nar_w[[stateName]]
            pNA * pNAR # TODO: alternatively, use external replication profile + cell$replicationModel$fractionS to set pNAR 
        })

        # solve the HMM
        hmm <- new_hmmEPTable(log(emissProbs), transProbs = log(transProbs), keys = windows$chrom) 
        cell$windows[[shapeKey]]$composite$stateI <<- keyedViterbi(hmm)

        # find and suppress mis-called segments in replicating cells that falsely nominate CNVs due to replication GC profile
        cell <<- suppressFalseCompositeCnvs(cell)
        cell$windows[[shapeKey]]$composite$HMM <<- cell$states[[shapeKey]][cell$windows[[shapeKey]]$composite$stateI, CN]
        cell$windows[[shapeKey]]$composite$NAR <<- cell$states[[shapeKey]][cell$windows[[shapeKey]]$composite$stateI, NAR]
    }
    iterateCompositeHmm(iteration = 1, transProb = 1e-10)
    tmp <- cell$windows$fitI
    cell$windows$fitI <- TRUE # resolve the HMM using the updated GC guidance from the first composite solution
    iterateCompositeHmm(iteration = 2, transProb = 1e-9)
    cell$windows$fitI <- tmp

    # go back and revise the sequential model to ensure that 
    if(is.null(cell$secondReplicationPass)){
        cell <- suppressFalseSequentialCnvs(cell)
        if(is.null(cell$windows[[shapeKey]]$sequential$NAR)){
            cell$replicationModel <- NULL
            cell$windows[[shapeKey]]$composite <- NULL
            name <- "masked"
            cell <- with(cell$windows[[shapeKey]], { 
                fitGcBias(cell, name, NR_wms, windowCN = sequential$HMM, col_w = defaultPointColor, composite = FALSE)
            })
            with(cell$windows[[shapeKey]], { plotWindows_cn(cell, name, FALSE, gc_w, NR_wms, TRUE, sequential$HMM) })
            cell$secondReplicationPass <- TRUE
            return(setReplicationModel(cell))
        }
    } else {
        cell$secondReplicationPass <- NULL
    }

    # finish up
    NR_wms <- cell$windows[[shapeKey]]$NR_wms
    hmm <- cell$windows[[shapeKey]]$composite$HMM
    plotWindows_cn(cell, "compositeByGc",  TRUE, gc_w, NR_wms, hmm = hmm, col = "GC")
    plotWindows_cn(cell, "compositeByCn",  TRUE, gc_w, NR_wms, hmm = hmm, col = "CN")
    plotWindows_cn(cell, "compositeByRep", TRUE, gc_w, NR_wms, hmm = hmm, col = "replication")
    finishShapeModel( setHMMProfiles(cell, composite = TRUE) )
}
#=====================================================================================

#=====================================================================================
# final post-processing and clean up of a cell object for tranmissal to app
#-------------------------------------------------------------------------------------
finishShapeModel <- function(cell){
    workingStep <<- "finishShapeModel"
    if(cell$replicationModel[[shapeKey]]$cellIsReplicating){
        cell$windows[[shapeKey]]$cnsd <- with(cell$windows[[shapeKey]], { min(sequential$cnsd, composite$cnsd) })
        cell$windows[[shapeKey]]$keep <- with(cell$windows[[shapeKey]], { sequential$keep || composite$keep })
        cell$windows[[shapeKey]]$fractionS <- cell$windows[[shapeKey]]$composite$fractionS
    } else {
        cell$windows[[shapeKey]]$cnsd <- cell$windows[[shapeKey]]$sequential$cnsd
        cell$windows[[shapeKey]]$keep <- cell$windows[[shapeKey]]$sequential$keep
        cell$windows[[shapeKey]]$fractionS <- 0.0
    }
    cell
}
finishCell <- function(cell){
    cell$cellIsReplicating <- cell$replicationModel[[shapeKey]]$cellIsReplicating
    cell$cnsd <- cell$windows[[shapeKey]]$cnsd
    cell$keep <- cell$windows[[shapeKey]]$keep
    cell$fractionS <- cell$windows[[shapeKey]]$fractionS
    for(key in c("nextWindows", "states", "tmp")) cell[[key]] <- NULL
    for(key in c("fitI")) cell$windows[[key]] <- NULL
    for(key in c("shape", "density")) {
        cell$windows$unshaped[[key]] <- NULL
        if(!is.null(cell$windows$shaped))  cell$windows$shaped[[key]]  <- NULL
        if(!is.null(cell$windows$batched)) cell$windows$batched[[key]] <- NULL
    }
    for(key in c("stateI")) {
        cell$windows$unshaped$sequential[[key]] <- NULL
        if(!is.null(cell$windows$shaped))  cell$windows$shaped$sequential[[key]]  <- NULL
        if(!is.null(cell$windows$batched)) cell$windows$batched$sequential[[key]] <- NULL
        if(cell$cellIsReplicating) {
            cell$windows$unshaped$composite[[key]] <- NULL
            if(!is.null(cell$windows$shaped))  cell$windows$shaped$composite[[key]] <- NULL
            if(!is.null(cell$windows$batched)) cell$windows$batched$composite[[key]] <- NULL
        }
    }
    for(key in c("NAME", "NAR", "FAR", "N_ALLELES")) {
        cell$segments$unshaped$sequential$CN[[key]] <- NULL
        if(!is.null(cell$segments$shaped))  cell$segments$shaped$sequential$CN[[key]] <- NULL
        if(!is.null(cell$segments$batched)) cell$segments$batched$sequential$CN[[key]] <- NULL
        if(cell$cellIsReplicating) {
            cell$segments$unshaped$composite$CN[[key]] <- NULL
            if(!is.null(cell$segments$shaped))  cell$segments$shaped$composite$CN[[key]] <- NULL
            if(!is.null(cell$segments$batched)) cell$segments$batched$composite$CN[[key]] <- NULL
        }
    }
    for(key in c("NAME", "CN",  "NAR", "N_ALLELES")) {
        cell$segments$unshaped$sequential$FAR[[key]] <- NULL
        if(!is.null(cell$segments$shaped))  cell$segments$shaped$sequential$FAR[[key]] <- NULL
        if(!is.null(cell$segments$batched)) cell$segments$batched$sequential$FAR[[key]] <- NULL
        if(cell$cellIsReplicating) {
            cell$segments$unshaped$composite$FAR[[key]] <- NULL
            if(!is.null(cell$segments$shaped))  cell$segments$shaped$composite$FAR[[key]] <- NULL
            if(!is.null(cell$segments$batched)) cell$segments$batched$composite$FAR[[key]] <- NULL
        }
    }
    cell    
}
#=====================================================================================

#=====================================================================================
# apply a further correction to replication squashed data that accounts for the shape
# of individual chromosome coverage profiles, which results from amplification artifacts
#-------------------------------------------------------------------------------------
getNR_wms <- function(NR_wm, shape) { # apply the chromosome shape correction to NR_wm
    if(shapeKey == "unshaped" || length(shape) == 1) return(NR_wm)
    NR_wms <- NR_wm / shape
    N_pre  <- sum(NR_wm,  na.rm = TRUE)
    N_post <- sum(NR_wms, na.rm = TRUE)
    NR_wms * N_pre / N_post # rescale to maintain total read count
}
processChromShape <- function(cell, windows, I){ # shape fit a single chromosome, checking that large CNVs aren't a result of shape effects
    fI <- cell$windows$fitI[I] # model is created from filtered windows
    nWorkingWindows <- sum(fI)
    if(nWorkingWindows < 10) return(rep(1, nrow(windows))) # insufficient points to fit
    i <- windows[, i]
    i_f <- i[fI] # collect values from this chrosomes fitted windows
    CN_wmshrf <- cell$tmp$CN_wmshr[I][fI]
    fitAllWindows <- function(){ # path used if all windows are either used or mask, no large CNVs to check
        fit <- loess(CN_wmshrf ~ i_f, span = 1, na.action = na.exclude)
        predict(fit, newdata = data.frame(i_f = i)) / cell$ploidy # final shape prediction occurs on ALL bins         
    }
    segments <- as.data.table(unclass(rle(cell$windows$unshaped$sequential$HMM[I][fI]))) # i.e., CN neutral and CNV spans
    if(nrow(segments) == 1) return(fitAllWindows()) # either euploid or completely aneuploid chromosome, fit as CN_wmshrf
    segments$hmm   <- segments$values
    segments$end   <- cumsum(segments$lengths)
    segments$start <- head(c(1, segments$end + 1), -1)
    chromModalCN   <- segments[, .(nWindows = sum(lengths)), by = "hmm"][order(-nWindows)][1, hmm]
    segments$modal <- segments$values == chromModalCN # the majority CN is used to fit the chromosome (not necessarily CN neutral)
    segments$mask  <- !segments$modal & segments$lengths < sum(segments$lengths) *  0.05
    segments$check <- !segments$modal & !segments$mask
    if(any(segments$mask)) segments[mask == TRUE, CN_wmshrf[start:end] <<- NA, by = start] # mask small non-modal, they will be ignored during shape fitting
    if(segments[, all(!check)]) return(fitAllWindows()) # chromosome only has small CNVs not included in fit, fit as CN_wmshrf_masked
    CN_wmsf_working <- copy(CN_wmshrf)
    CN_wmsrf <- cell$tmp$CN_wmsr[I][fI]
    segments <- segments[check == TRUE][order(-lengths)] # check the largest segmental CNVs first
    segments[, CN_wmsf_working[start:end] <<- NA, by = start] # initially, mask all segments we will check (filled in again below)
    checkCnvFit <- function(){
        fit <- loess(CN_wmsf_working ~ i_f, span = 1, na.action = na.exclude)
        delta_norm <- resid(fit) / predict(fit) # must scale the residuals to account for CN spreading
        sqrt(sum(delta_norm ** 2, na.rm = TRUE)) # use RMSD to compare +/- CNV fits 
    }
    for(segmentI in 1:nrow(segments)){
        segments[segmentI, CN_wmsf_working[start:end] <<- CN_wmshrf[start:end]]
        rmsd_keepCnv <- checkCnvFit() # this model DOES correct for the CNV, thus assumes it is real
        segments[segmentI, CN_wmsf_working[start:end] <<- CN_wmsrf[start:end]]
        rmsd_rejectCnv <- checkCnvFit() # this model DOES NOT correct for the CNV, instead calculates the shape of the CN-shifted windows
        if(rmsd_keepCnv < rmsd_rejectCnv){ # use the model that has the best loess fit with least residuals
            segments[segmentI, CN_wmsf_working[start:end] <<- CN_wmshrf[start:end]] # this CNV span is alway filled before checking the next one
        }
    }
    fit <- loess(CN_wmsf_working ~ i_f, span = 1, na.action = na.exclude) # return the shape after making decision on all CNV spans
    predict(fit, newdata = data.frame(i_f = i)) / cell$ploidy
}
setChromShapes <- function(cell){
    workingStep <<- "setChromShapes"
    windows <- windows[[cell$windowPower + 1]]
    gc_w <- windows[, gc_fraction]    
    fitI <- cell$windows$fitI
    NA_ <- with(cell$windows$unshaped$sequential, { HMM + NAR })
    cell$tmp$NR_wmshr <- with(cell, { windows$unshaped$NR_wms * ploidy / (windows$unshaped$sequential$HMM + windows$unshaped$sequential$NAR) }) # thus, correct counts toward ploidy
    cell$tmp$NR_wmsr  <- with(cell, { windows$unshaped$NR_wms * windows$unshaped$sequential$HMM / NA_ }) # correct toward unreplicated, CNVs persist
    RPA <- predict(cell$replicationModel$unshaped$gc_fit, gc_w) # this is an unsquashed GC model corrected for both HMM and NAR
    cell$tmp$CN_wmshr <- cell$tmp$NR_wmshr / RPA
    cell$tmp$CN_wmsr  <- cell$tmp$NR_wmsr  / RPA
    shape <- windows[, { # shape chromosomes one at a time
        processChromShape(cell, .SD, .I)
    }, by = "chrom"][[2]]
    plotWindows_cn(cell, "shapes", FALSE, gc_w, cell$tmp$NR_wmshr, windowI = TRUE, 
                   gc_fit = cell$replicationModel$unshaped$gc_fit, shape = shape,
                   col = defaultPointColor)
    cell$nextWindows <- NULL # go back to the beginning to solve GC and replication, now with the benefit of chromosome shapes
    cell$tmp <- NULL
    shapeKey <<- "shaped" 
    cell$windows$shaped <- fitCellWindows(cell$windowPower, windows, cell$windows$unshaped$NR_wms, shape)    
    cell
}
#=====================================================================================

#=====================================================================================
# apply a further correction to normalize for batch effects accross a given sample of cells
#-------------------------------------------------------------------------------------
setBatchShapes <- function(cell){
    workingStep <<- "setBatchShapes"
    windows <- windows[[cell$windowPower + 1]]
    cell$nextWindows <- NULL # go back to the beginning to solve GC and replication, now with the benefit of chromosome shapes
    cell$tmp <- NULL
    workingShapeKey <- shapeKey
    shapeKey <<- "batched"     
    cell$windows$batched <- fitCellWindows(cell$windowPower, windows, cell$windows[[workingShapeKey]]$NR_wms, cell$windows$batched$shape)
    plotCountDensity(cell)
    plotNumber <<- plotNumber + 1
    cell 
}
#=====================================================================================

#=====================================================================================
# main target function that launches the process
#-------------------------------------------------------------------------------------
fitCell_1 <- function(cell_id){ # creates the unshaped and shaped, i.e., cell-specific, fits
    plotNumber <<- 1 # so developer plots are ordered in lists as they were created
    workingStep <<- "fitCell_1"
    shapeKey <<- "unshaped"
    cell <- list(
        cell_id  = cell_id,
        badCell  = FALSE,
        modal_NA = env$PLOIDY, # subject to change if late S-phase
        ploidy   = env$PLOIDY,  # never changed by this script, always env$PLOIDY
        replicationModel = list(),
        states = list(),
        segments = list()
    )
    tryCatch({    
        cell$minWindowCount <- getMinWindowCount(env$PLOIDY)
        cell$minWindowPower <- getMinWindowPower(cell_id, cell$minWindowCount)
        cell <- setCellWindows(cell)
        if(cell$badCell){
            cell$windows$unshaped$density <- NULL
        } else {
            if(env$SHAPE_CORRECTION %in% c('cell', 'both')){
                cell <- setChromShapes(cell)
                cell <- setChromCN(cell)
            }
            cell <- finishCell(cell)            
        }
    }, error = function(e) {
        message(paste("error at fitCell_1 step:", workingStep))
        message(paste("cell_id", cell_id))
        print(e)
        cell$badCell <- TRUE
    })
    cell
}
fitCell_2 <- function(cell){ # creates the batched fits per cell, after running normalizeBatch on a set of cells from a sample
    plotNumber <<- 100
    workingStep <<- "fitCell_2"    
    if(cell$badCell || is.null(cell$windows$batched)) return(cell)
    tryCatch({    
        shapeKey <<- if(env$SHAPE_CORRECTION %in% c('cell', 'both')) "shaped" else "unshaped"
        cell <- setBatchShapes(cell)
        cell <- setChromCN(cell)
        plotBatchEffect(cell)
        cell <- finishCell(cell)
    }, error = function(e) {
        message(paste("error at fitCell_2 step:", workingStep))
        message(paste("cell_id", cell$cell_id))
        print(e)
        cell$badCell <- TRUE
    })
    cell
}
#=====================================================================================
