#=====================================================================================
# throughout this script and pipeline, the following naming conventions are used
#-------------------------------------------------------------------------------------
#   variable names:
#     NR  = Number of Reads, i.e., a read count (can be fractional due to fragment splitting across bins)
#     ER  = Expected Reads, i.e., the nominal read count that corresponds to a specific copy number and/or state
#     RPA = Reads Per Allele, i.e., the count expected for CN == 1, i.e., ER_ploidy / ploidy
#     CN  = Copy Number, PRIOR to any replication, can be fractional estimates as NR / ER
#     HMM = Hidden Markov Model, i.e., quantal CN values determined from a model of NR + ER
#     NA_ = Number of Alleles, i.e., how many DNA copies are present considering BOTH CN and replication (underscore since NA is reserved)
#     NAR = Number of Alleles Replicated, thus, could be c(0,1,2) for CN==2, to yield NA=c(2,3,4)
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
#     h = counts that have been rescaled to an HMM, e.g., NR * ploidy / HMM
#     r = counts that have been rescaled to predicted replication states
#     s = counts that have been rescaled based on chromosome shape biases (see below)
#-------------------------------------------------------------------------------------
# additional considerations:
#   - "ploidy" is the expected autosome copy number prior to replication, either 1 or 2 at present
#   - "modal" and "peak" are largely interchangeable and refer to the mode of a set of NR_w values
#   - "windowSize" is the number of 20kb bins that comprise a cell's windows, can vary between cells
#   - "windowPower" = log2(windowSize), i.e., allowed window sizes increment by doubling
#   - "GC bias" is the dependence of NR values on window GC base content
#   - "squashed" refers to a GC bias correction that does not consider replication, 
#                i.e., that masks replication-based differences in NR as a result of GC bias correction
#   = "shape" is an artifact where a chromosome's counts are not flat but tilted, S-shaped, etc.
#                 shape effects vary considerably with the amplification method
#   - "fractionS" is the total fraction of alleles across the genome that have been replicated (FAR), 
#                 e.g., at fractionS==0.5, half of the genome's bases (not windows!) have been replicated
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
getAgg <- function(NR){ # mostly used for debugging
    x <- round(NR, 0)
    agg <- aggregate(x, list(x), length)
    names(agg) <- c("NR", "N")
    agg$freq <- agg$N / sum(agg$N)
    agg
}
#=====================================================================================

#=====================================================================================
# plotting functions
#-------------------------------------------------------------------------------------
plotNumber <- 1
defaultPointColor <- rgb(0, 0, 0, 0.1)
plotsDir <- file.path(env$PLOTS_DIR, 'scCNV')
if(!dir.exists(plotsDir)) dir.create(plotsDir)
# plotsDir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV/plots"
saveCellPlot <- function(cell, name, fn, width = 2, height = 2){ # save a plot for future assembly and app
    plotDir <- file.path(plotsDir, cell$cell_id)
    if(!dir.exists(plotDir)) dir.create(plotDir)
    filename <- paste(cell$cell_id, letters[plotNumber], name, "png", sep = ".")
    pngFile <- file.path(plotDir, filename)
    png(pngFile, width = width, height = height, units = "in", pointsize = 7, res = 300, type = "cairo")
    fn()
    dev.off()
    plotNumber <<- plotNumber + 1
}
getRepColor <- function(cell, default = NULL){ # after replication fitting, color points by replication state
    if(is.null(cell$replicationModel)) {
        if(is.null(default)) defaultPointColor else default
    } else {
        cols <- c(rgb(0, 0, 1, 0.15), rgb(0, 1, 0, 0.15), rgb(1, 0, 0, 0.15), rgb(0, 0.75, 0.75, 0.15), rgb(1, 0, 1, 0.15))
        col <- cols[cell$windows$NAR + 1]
        col[is.na(col)] <- defaultPointColor
        col
    }
}
getCnColor <- function(cell){
    cols <- c(defaultPointColor, rgb(0, 0, 1, 0.15), rgb(0, 1, 0, 0.15), rgb(1, 0, 0, 0.15), rgb(0, 0.75, 0.75, 0.15), rgb(1, 0, 1, 0.15)) 
    col <- cols[cell$states[cell$windows$stateI, CN + 1]]
    col[is.na(col)] <- defaultPointColor
    col
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
plotGcBias <- function(cell, name, gc_w, NR_wm, windowCN = NULL){ # NR vs. fraction GC, sometimes with adjustments for CN and/or replication
    saveCellPlot(cell, paste("gc", name, sep = "_"), function(){
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        plot(cell$gc_fit, gc_w, NR_wm, cell$ploidy, col = getRepColor(cell, NA), binCN = windowCN)
    })
}
plotWindows_counts <- function(cell, name, NR_wm, windowI = TRUE, col = defaultPointColor, shape = NULL){ # plot window counts by chromosome coordinates
    saveCellPlot(cell, paste("NR_wm", name, sep = "_"), function(){
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        coord <- getWindowCoordinates(cell, windowI)
        if(is.character(col) && col[1] == "useRepColor") col <- getRepColor(cell)
        if(is.character(col) && col[1] == "useCnColor")  col <- getCnColor(cell)
        plot(NA, NA, xlim = coord$xlim, ylim = c(0, cell$windows$ER_modal_NA * 3),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "# of Reads")
        addChromLabels(coord)
        if(!is.null(cell$windows$RPA)) abline(h = cell$windows$RPA * 0:10, col = "grey")        
        points(coord$x, NR_wm, pch = 19, cex = 0.3, col = col)
        if(!is.null(shape)) points(coord$x, shape, pch = 16, cex = 0.5, col = "blue")
    }, width = 6)
}
plotWindows_cn <- function(cell, name, gc_w, NR_wm, windowI = TRUE, hmm = NULL){ # plot window CN by chromosome coordinates
    saveCellPlot(cell, paste("cn", name, sep = "_"), function(){
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        RPA <- predict(cell$gc_fit, gc_w, type = 'adjustedPeak')
        cn <- NR_wm / RPA
        coord <- getWindowCoordinates(cell, windowI)
        col <- getRepColor(cell)
        plot(NA, NA, xlim = coord$xlim, ylim = c(0, maxModelCN + 1),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "Copy Number")
        addChromLabels(coord)
        abline(h = 0:(maxModelCN + 1), col = "grey")
        points(coord$x, cn, pch = 19, cex = 0.3, col = col)
        if(!is.null(hmm)) lines(coord$x, hmm, col = "red3")
    }, width = 6)
}
plotWindows_FAR <- function(cell, name){ # plot window CN by chromosome coordinates
    saveCellPlot(cell, paste("FAR", name, sep = "_"), function(){
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        coord <- getWindowCoordinates(cell, TRUE)
        col <- getRepColor(cell)
        plot(NA, NA, xlim = coord$xlim, ylim = c(-0.5, 2),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "Fraction Replicated")
        addChromLabels(coord)
        abline(h = 0:1, col = "grey")
        points(coord$x, cell$windows$FAR, pch = 19, cex = 0.3, col = col)
        # lines(coord$x, cell$windows$NAR / cell$windows$HMM, col = "red3") # too busy...
    }, width = 6)
}
plotCellFit <- function(cell, name, model, bestModelByType){ # distributions of read counts with best-fitted replication models
    saveCellPlot(cell, paste("repFit", name, sep = "_"), function(){
        NR_a <- round(cell$windows$NR_wmhfl, 0) # here, a is for "actual"
        N <- length(NR_a)
        getRandomWindows <- function(model){ with(model, {
            weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, fractionS) # see utilities/create_fractionS_table.R
            N_R0 <- round(N * weights[1], 0) # no  alleles replicated
            N_R1 <- round(N * weights[2], 0) # one allele  replicated
            N_R2 <- round(N * weights[3], 0)   # two alleles replicated
            round(c(
                rnbinom(N_R0, mu = ER_unreplicated,       size = theta),
                rnbinom(N_R1, mu = ER_unreplicated * 1.5, size = theta),
                rnbinom(N_R2, mu = ER_unreplicated * 2,   size = theta)
            ), 0)              
        })}
        NR_m <- getRandomWindows(model)
        NR_ad <- density(NR_a)
        NR_md <- density(NR_m)
        plot(NR_ad, typ="l", col = "black", lwd = 2,
             xlim = c(0, max(NR_ad$x, NR_md$x)), ylim = c(0, max(NR_ad$y, NR_md$y)))
        abline(v = cell$windows$ER_modal_NA * c(1/2, 1, 2))
        bestModelByType[, {
            NR_x <- getRandomWindows(.SD)
            col <- modelTypeColors[[modelType]]
            lwd <- if(modelType == model$modelType) 2 else 1
            lines(density(NR_x), col = col, lwd = lwd)
            abline(v = ER_unreplicated * 1:2, lty = 2, col = col)
        }, by = "modelType"]
    }, width = 3, height = 3)
}
#=====================================================================================

#=====================================================================================
# functions that help adjust the window size to provide reliable copy number state shifts
#-------------------------------------------------------------------------------------
# determine the least number of reads per window for a cell's data to support a robust HMM
# places ~96% of windows within (modal_NA +/- 0.5) * replicationFactor==[1,2]
# sapply(1:4, getMinWindowCount) => 16 64 144 256; thus, usually 64 reads/window for interphase diploid cells
getMinWindowCount <- function(modal_NA){
    ploidyFactor <- (modal_NA + 0.5) / modal_NA
    (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# require that most windows have non-zero count, with allowance for rare homozygous losses and chrY
# similarly, the window median must have a usable count
# reject a windowPower for a cell that cannot achieve these metrics
checkForExcessZeros <- function(cell_id, windowPower){
    windowSize <- 2 ** windowPower
    NR_w <- collapseVector(raw_counts[[cell_id]], windowSize)
    NR_wf <- NR_w[!is.na(NR_w)]
    mean(NR_wf < 1) > 0.2 || # return TRUE (reject) if too many zero-count windows
    median(NR_wf) < 5        # or if median is too low such that most windows have too few counts
}

# set the per-cell minimum window size as the number of bins needed to obtain a mean raw count >=minWindowCount
# cells with too-low counts will still be analyzed at MAX_WINDOW_POWER but will likely later fail QC
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
    windowPower
}

# use Q95 of window-to-window count deltas (NOT of raw window values) to determine 
# the window size that puts the majority of window counts between CN = ploidy +/- 1
# for most good cells, this default windowPower is > minWindowPower due to overdispersion
normLagDiffQ_quantile  <- 0.5   # determined empirically
normLagDiffQ_threshold <- 0.175 # determined empirically
getNormLagDiffQ <- function(NR){ # 95%ile of the magnitude of the window-to-window count difference, normalized to local mean count
    lagDiff <- abs(diff(NR)) 
    lagMean <- (head(NR, -1) + tail(NR, -1)) / 2
    quantile(lagDiff / lagMean, normLagDiffQ_quantile, na.rm = TRUE)
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
    NR_wmf <- excludeOutliers(NR_wm, min = 1) # thus, chrom-to-chrom, baseline-to-CNV, and rare bad windows won't have undue influence
    ER_modal_NA <- peakValue(NR_wmf) # ER_modal_NA may NOT be the same as ER_ploidy 
    list( # one of these is chosen to become the working cell object
        NR_wm       = NR_wm,        
        ER_modal_NA = ER_modal_NA,
        ER_ploidy   = ER_modal_NA,
        RPA         = ER_modal_NA / env$PLOIDY, # later we may learn RPA is 2-fold off for replicating cells,
        normLagDiffQ = getNormLagDiffQ(NR_wmf)
    )
}
setCellWindows <- function(cell){ # establish the optimal window power for a cell
    # message("setCellWindows")
    cell$windowPower <- cell$minWindowPower
    cell$windows <- getCellWindows(cell)
    while(cell$windows$normLagDiffQ > normLagDiffQ_threshold && # stop when we get a nice, appropriately tight distribution
          cell$windowPower < env$MAX_WINDOW_POWER){
        cell$windowPower <- cell$windowPower + 1
        cell$windows <- getCellWindows(cell)
    }
    setChromCN(cell)
}
#=====================================================================================

#=====================================================================================
# shared cell fitting functions, for GC bias effects and copy number by HMM
#-------------------------------------------------------------------------------------
fitGcBias <- function(cell, name, NR_wm, windowI = TRUE, chroms = NULL, windowCN = NULL, col_w = NULL){
    windows <- windows[[cell$windowPower + 1]]
    if(!is.null(chroms)) windowI <- windowI & windows[, chrom %in% chroms]        
    gc_wf <- windows[windowI, gc_fraction]
    NR_wmf <- NR_wm[windowI]
    if(is.null(windowCN)) windowCN <- cell$ploidy
    if(length(windowCN) > 1 && length(windowCN) != length(NR_wmf)) windowCN <- windowCN[windowI]
    cell$gc_fit <- new_nbinomCountsGC2(NR_wmf, gc_wf, binCN = windowCN)
    plotWindows_counts(cell, name, NR_wm, col = col_w)    
    plotGcBias(cell, name, gc_wf, NR_wmf, windowCN = windowCN)
    cell
}
solveSimpleHMM <- function(cell, name, NR_wm, windowI = TRUE){
    windows <- windows[[cell$windowPower + 1]]
    gc_wf <- windows[windowI, gc_fraction]
    NR_wmf <- NR_wm[windowI]
    cell$windows$HMM <- viterbi(
        cell$gc_fit, 
        NR_wmf, # note, could simply be NR_wm if windowI == TRUE
        gc_wf, 
        maxCN = maxModelCN, 
        chroms = windows[windowI == TRUE, chrom], 
        asRle = FALSE
    )$cn
    plotWindows_cn(cell, name, gc_wf, NR_wmf, windowI, cell$windows$HMM)
    cell
}
#=====================================================================================

#=====================================================================================
# use a replication-unaware algorithm that squashes the replication GC bias effect
# to obtain an initial copy number estimates across the genome
#-------------------------------------------------------------------------------------
maxModelCN <- 5
CNs <- 0:maxModelCN # CN==5 is not trustworthy, the true value could be (much) higher
setChromCN <- function(cell){ 
    # message("setChromCN")
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # determine a first crude aneuploidy copy number estimate for each whole chromosome
    NR_wmf <- excludeOutliers(cell$windows$NR_wm, min = 2)
    if(length(NR_wmf) > 1e4) NR_wmf <- sample(NR_wmf, 1e4)
    RPA_wmf <- NR_wmf / cell$ploidy
    NR_wmf_cn <- lapply(CNs, "*", RPA_wmf)
    cell$windows$CN_c <- lapply(chroms, function(chrom){
        NR_wmcl <- excludeOutliers(cell$windows$NR_wm[windows$chrom == chrom], min = 2)
        if(length(NR_wmcl) <= 1 || median(NR_wmcl) <= 1) return(0)
        p <- sapply(NR_wmf_cn, function(NR_g) wilcox.test(NR_g, NR_wmcl, exact = FALSE)$p.value)
        CNs[which.max(p)] # the copy number that gives the best distribution match between chrom and genome
    })
    names(cell$windows$CN_c) <- chroms

    # use the crude estimates to select chromosomes for performing a proper GC bias fit and CN HMM
    chromFitI <- unlist(cell$windows$CN_c) == cell$ploidy
    col_w <- sapply(windows$chrom, function(chrom) {
        if(cell$windows$CN_c[[chrom]] == cell$ploidy) defaultPointColor else rgb(1, 0, 0, 0.1)
    })
    cell <- fitGcBias(cell, "input", cell$windows$NR_wm, chroms = chroms[chromFitI], col_w = col_w)
    cell <- solveSimpleHMM(cell, "squashed", cell$windows$NR_wm)
    setReplicationModel(cell) # pass the baton
}
#=====================================================================================

#=====================================================================================
# establish the common parameters of a set of replication models to be applied to individual cells
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
getModelLogLikelihood <- function(P_R0, P_R1, P_R2, fractionS){ # P_X from density functions, to be weighted by fractionS
    weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, fractionS) # see utilities/create_fractionS_table.R
    x <- log(
        P_R0 * weights[1] + # zero alleles replicated
        P_R1 * weights[2] + # one  allele  replicated
        P_R2 * weights[3]   # two  alleles replicated
    )
    x[x == -Inf] <- -800 # suppress outliers      
    sum(x)
}
getModelER_unreplicated <- function(modelType, cell){ # set read count expectations (ER) for a combination of modelType and cell$ER_modal_NA
    with(cell$windows, { switch(
        modelType,
        notReplicating     = ER_modal_NA,
        peakIsUnreplicated = ER_modal_NA,    # the modal window count corresponds to unreplicated DNA
        peakIsReplicated   = ER_modal_NA / 2 # the modal window count corresponds to   replicated DNA
    )})
}
fillModelLikelihoods <- function(modelType, fractionS, cell){

    # use a rich subset of the data to fit the candidate model to the data as best as possible
    # models are defined by replication peak type and fractionS increments
    ER_unreplicated <- getModelER_unreplicated(modelType, cell)
    optFn <- function(par){ 
        LL <- with(cell$windows, { -getModelLogLikelihood(
            dnbinom(NR_wmhfld, mu = par[1],       size = par[2]), 
            dnbinom(NR_wmhfld, mu = par[1] * 1.5, size = par[2]),
            dnbinom(NR_wmhfld, mu = par[1] * 2,   size = par[2]),
            fractionS
        )})
        LL
    }
    par <- optim( # by fitting both ER_unreplicated and theta, we arrive at mu and size parameters for nbinom
        c(ER_unreplicated = ER_unreplicated, theta = 7), # adjusted from our initial estimate made from the peak, i.e., mode
        optFn,
        method = "L-BFGS-B",
        lower = c(ER_unreplicated = ER_unreplicated * 0.9, theta = 1),
        upper = c(ER_unreplicated = ER_unreplicated * 1.1, theta = 1000)
    )$par

    # use those values to establish the log likelihood for this modelType + fractionS over all peaks
    par[3] <- with(cell$windows, { getModelLogLikelihood(
        dnbinom(NR_wmhfl, mu = par[1],       size = par[2]), 
        dnbinom(NR_wmhfl, mu = par[1] * 1.5, size = par[2]), 
        dnbinom(NR_wmhfl, mu = par[1] * 2,   size = par[2]),
        fractionS
    ) })
    par
}
setBestReplicationModel <- function(cell, name){

    # fit all required models
    models <- copy(modelValues)
    if(!is.null(cell$replicationModel)){ # when refitting for ER_modal_NA and theta, just use the single established model
        models <- with(cell, { models[ modelType == replicationModel$modelType & fractionS == replicationModel$fractionS ] })
    }
    models <- cbind(models, t(models[, mapply(function(modelType, fractionS){
        fillModelLikelihoods(modelType, fractionS, cell)
    }, modelType, fractionS)]))
    setnames(models, c("modelType", "fractionS", "ER_unreplicated", "theta", "logLikelihood"))

    # select the best model
    bestModelByType <- models[, {
        i <- which.max(logLikelihood)
        .SD[i]
    }, by = "modelType"]
    cell$replicationModel <- as.list(models[which.max(logLikelihood)])
    plotCellFit(cell, name, cell$replicationModel, bestModelByType) # based on NR_wmhfl

    # calculate derived model parameters
    peakIsReplicated <- cell$replicationModel$modelType == "peakIsReplicated"
    cell$modal_NA  <- cell$ploidy * (if(peakIsReplicated) 2 else 1)
    cell$windows$ER_ploidy <- cell$windows$ER_modal_NA / (if(peakIsReplicated) 2 else 1)    
    cell$windows$RPA <- cell$windows$ER_modal_NA / cell$modal_NA 
    cell$cellIsReplicating <- cell$replicationModel$modelType != "notReplicating"
    cell
}
getRepEmissProbs <- function(cell){
    model <- cell$replicationModel
    weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, model$fractionS)
    with(cell$windows, { log(as.matrix(data.table( # here, model$ER_unreplicated has been fitted to reflect mu, not peak/mode
        dnbinom(NR_wmhf, mu = model$ER_unreplicated,       size = model$theta) * weights[1], 
        dnbinom(NR_wmhf, mu = model$ER_unreplicated * 1.5, size = model$theta) * weights[2], 
        dnbinom(NR_wmhf, mu = model$ER_unreplicated * 2,   size = model$theta) * weights[3]
    ))) })
}
#=====================================================================================

#=====================================================================================
# fit a cell's data to establish the initial parameters of it's replication model
#-------------------------------------------------------------------------------------
setReplicationModel <- function(cell){
    # message("setReplicationModel")
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # use CNV-corrected counts to establish cell-specific likelihoods for replication models
    cell$windows$NR_wmh <- with(cell, { windows$NR_wm * ploidy / windows$HMM }) 
    repFitI <- with(cell$windows, { HMM > 0 & HMM < maxModelCN & !is.na(NR_wmh) }) 
    cell$windows$NR_wmhf <- round(cell$windows$NR_wmh[repFitI], 0)
    cell$windows$NR_wmhfl <- excludeOutliers(cell$windows$NR_wmhf)
    cell$windows$NR_wmhfld <- with(cell$windows, { if(length(NR_wmhfl) > 1e4) sample(NR_wmhfl, 1e4) else NR_wmhfl }) 

    # optimize theta for three distinct types of models, at multiple fixed fractionS values
    # and ER_(un)replicated values determine by model type
    # pick the best model type, with its best parameters
    cell <- setBestReplicationModel(cell, "squashed")

    # solve an HMM to establish the (un)replicated genome spans for adjusting window NR prior to GC fit
    # message("solve rep HMM")
    emissProbs <- getRepEmissProbs(cell)
    hmm <- new_hmmEPTable(emissProbs, transProb = 1e-1, keys = windows$chrom[repFitI])
    cell$windows$NAR <- {
        NAR <- keyedViterbi(hmm) - 1 # based on NR_wmhf
        round(NAR * cell$windows$HMM[repFitI] / cell$ploidy, 0) # rescale NAR to hmm
    }
    plotWindows_counts(cell, "repHMM", cell$windows$NR_wm[repFitI], windowI = repFitI, col = "useRepColor") 
    setChromShapes(cell, repFitI) # pass the baton
}
#=====================================================================================

#=====================================================================================
# apply a further correction to replication squashed data that accounts for the shape
# of individual chromosome coverage profiles, which results from amplification artifacts
#-------------------------------------------------------------------------------------
setChromShapes <- function(cell, shapeFitI){ # shapeFitI is the same as repFitI
    # message("setChromShapes")
    windows <- windows[[cell$windowPower + 1]]

    # fit a quadratic to each chromosome and use to normalize counts, i.e., stabilize the NR baseline
    NAR <- rep(NA, nrow(windows)) # begin with a list of all windows, only some of which have NAR data
    NAR[shapeFitI] <- cell$windows$NAR
    cell$windows$NR_wmhr <- with(cell, { windows$NR_wm * modal_NA / (windows$HMM + NAR) }) # thus, correct counts toward ploidy
    cell$windows$shape <- windows[, { # shape is relative to ER_modal_NA
        sfI  <- shapeFitI[.I] # model is created from filtered bins
        if(sum(sfI) < 10) rep(1, .N) else { # insufficient points to fit
            i_f  <- i[sfI]
            i2_f <- i2[sfI]
            NR_wmhrf <- cell$windows$NR_wmhr[.I][sfI]
            if(all(is.na(NR_wmhrf))) rep(1, .N)
            fit <- lm(NR_wmhrf ~ i_f + i2_f)
            predict(fit, newdata = data.frame(i_f = i, i2_f = i2)) # final shape prediction occurs on ALL bins 
        }
    }, by = "chrom"][[2]]
    cell$windows$shapeCorrection <- with(cell$windows, { shape / ER_modal_NA })     
    cell$windows$NR_wms <- with(cell$windows, { NR_wm / shapeCorrection })
    cell$windows$NR_wms <- with(cell$windows, { NR_wms * sum(NR_wm, na.rm = TRUE) / sum(NR_wms, na.rm = TRUE) }) # maintain total cell weights
    plotWindows_counts(cell, "shapes", cell$windows$NR_wmhr[shapeFitI],
                       windowI = shapeFitI, shape = cell$windows$shape[shapeFitI])

    # re-solve the GC bias fit using the reshaped NR values
    cell$windows$NR_wmshr <- with(cell, { windows$NR_wms * ploidy / (windows$HMM + NAR) }) # thus, correct counts toward ploidy
    col_w <- ifelse(shapeFitI, defaultPointColor, rgb(1, 0, 0, 0.1))
    cell <- fitGcBias(cell, "reshaped", cell$windows$NR_wmshr, windowI = shapeFitI, col_w = col_w)
    plotWindows_counts(cell, "reshaped", cell$windows$NR_wms[shapeFitI], windowI = shapeFitI, col = "useRepColor")  
    solveCompositeHmm(cell, shapeFitI) # pass the baton
}
#=====================================================================================

#=====================================================================================
# support functions for replication state profiles as a function of fractionS + window GC
# i.e., that add a (P_NAR | GC) weight to the composite HMM for replicating cells
# derived from the intial fit of each individual cell (TODO: incorporate external replication weights?)
#-------------------------------------------------------------------------------------
nGcSteps <- 50
gcIndices  <- 1:nGcSteps
P_rep_fs_gc <- list() # for collecting replicating cell profiles for later aggregation
getP_replicated_wf <- function(cell, gc_wf) with(cell, {
    gci_wf <- round(gc_wf * nGcSteps, 0) # vector of filtered window GC indices
    sapply(gcIndices, function(gci){ # all bins have at least somelikelihood of being called replicating
        NAR <- windows$NAR[gci_wf == gci] # aggregate all windows with a given GC index
        if(length(NAR) == 0) return(NA)
        agg <- aggregate(NAR, list(NAR), length)
        getN <- function(nar){
            n <- agg[agg[[1]] == nar, 2]
            if(length(n) == 0) 0 else n
        }
        pNAR <- sapply(0:ploidy, getN) / sum(agg[[2]]) # NAR's as observed for ploidy values in initial modeling
        minPNar <- 0.01
        if(ploidy == 1){
            sapply(0:4, function(cn) pmax(minPNar,
                     if(cn == 0) 1    # one state, replication not meaningful
                else if(cn == 1) pNAR # CN==1, same as the NAR model
                else if(cn == 2) c(pNAR[1], minPNar, pNAR[2]) # CN=2,3,4, interpolate intermediate values
                else if(cn == 3) c(pNAR[1], minPNar, minPNar, pNAR[2])
                else             c(pNAR[1], minPNar, minPNar, minPNar, pNAR[2])
            ))
        } else { # ploidy == 2, typical
            sapply(0:4, function(cn) pmax(minPNar,
                     if(cn == 0) 1                       # one state, replication not meaningful
                else if(cn == 1) c(pNAR[1], 1 - pNAR[1]) # two states, unreplicated or replicated
                else if(cn == 2) pNAR                    # CN==2, same as the NAR model
                else if(cn == 3) c(pNAR[1], pNAR[2], pNAR[2], pNAR[3]) # CN=3,4
                else             c(pNAR[1], pNAR[2], pNAR[2], pNAR[2], pNAR[3])
            ))
        }
    })
})
getP_replicated_w <- function(cell, gc_w)  with(cell$replicationModel, {
    gci_w <- round(gc_w * nGcSteps, 0) # vector of all window GC indices
    is <- 1:length(gci_w)    
    nonNaGci <- which(!is.na(P_replicated_wf)) # GC indices with sufficient filtered windows to estimate P_unreplicated_wf
    sapply(is, function(i){ # return the appropriate GCI per window (may have been adjusted from the actual GCI)
        gci <- gci_w[i]
        x <- P_replicated_wf[gci]
        if(length(x) > 0 && !is.na(x)) return(gci)
        delta <- abs(nonNaGci - gci)
        nonNaGci[which(delta == min(delta))[1]]
    })
})
commitReplicationProfile <- function(cell, gc_wf) with(cell, { # keep track of GC replication bias across cells
    if(!cellIsReplicating) return(NULL)
    gci_wf <- round(gc_wf * nGcSteps, 0) 
    P_rep_fs_gc[[cell_id]] <<- list(
        cell_id = cell_id,
        windowPower = windowPower,
        modelType = replicationModel$modelType,
        fractionS = replicationModel$fractionS,
        profile = sapply(gcIndices, function(gci){
            NAR <- replicationModel$NAR[gci_wf == gci]
            if(length(NAR) == 0) return(rep(NA, ploidy + 2))
            agg <- aggregate(NAR, list(NAR), length)
            N <- sum(agg[[2]])
            c(N, sapply(0:ploidy, function(nar) {
                n <- agg[agg[[1]] == nar, 2]
                if(length(n) == 0) 0 else n
            }) / N)
        })
    )  
})
# plotReplicationProfiles <- function(){
#     values <- c("N", "P_NAR0", "P_NAR1", "P_NAR2")
#     fractionS_int <- as.integer(seq(0.05, 0.95, 0.025) * 1000)
#     colfunc <- colorRampPalette(c("blue", "red"))
#     colors <- colfunc(length(fractionS_int))
#     weights <- sapply(P_rep_fs_gc, function(cell) cell$profile[1,])
#     for(j in 1:length(values)){ # one plot per value type
#         message(values[j])
#         saveDevPlot(paste("P_rep_fs_gc", values[j], sep = "."), function(){
#             x  <- gcIndices
#             ys <- sapply(P_rep_fs_gc, function(cell) cell$profile[j,])
#             plot(
#                 NA, NA, 
#                 xlim = c(0.3, 0.7), ylim = range(ys, na.rm = TRUE), 
#                 xlab = "Fraction GC", ylab = values[j]
#             ) 
#             for(k in 1:ncol(ys)) { # once trace per cell
#                 fractionS <- 
#                 P_rep_fs_gc[[k]]$replicationModel$fractionS
#                 col <- colors[fractionS_int == as.integer(fractionS * 1000)]
#                 y <- ys[, k]
#                 lines(x / nGcSteps, y, col = col)  
#                 points(x / nGcSteps, y, col = col)               
#             }
#         })
#     }
# }
#=====================================================================================

#=====================================================================================
# use the replication model and reshaped GC bias to solve a multi-state HMM with CN and replication components
#   use cell$NR_wms to normalize to chromosome shape, i.e., amplification skew
#   use cell$gc_fit to determine ER per state as gc_fit(GC) * (CN + NAR = NA)
#   use P_unreplicated_wf to weight windows for timed replication potential based on GC content
# none of these are deterministic for a window, they are only trends determined from aggregated data
# thus, the model output is not dictated or forced at any window by upstream preparative work
#-------------------------------------------------------------------------------------
compositeStates <- fread(file.path(env$ACTION_DIR, "composite_states.csv"))
solveCompositeHmm <- function(cell, shapeFitI){ # shapeFitI same as repFitI
    # message("solveCompositeHmm")
    windows <- windows[[cell$windowPower + 1]]

    # establish the cell specific values for (P_unreplicated | GC)
    gc_w <- windows[, gc_fraction]
    gc_wf <- gc_w[shapeFitI]
    cell$replicationModel$P_replicated_wf <- getP_replicated_wf(cell, gc_wf)
    P_replicated_wf <- cell$replicationModel$P_replicated_wf # since cell is reset below
    P_replicated_w <- getP_replicated_w(cell, gc_w) # a vector of reference GCI values

    # re-fit ER_modal_NA (and thus RPA) and theta for our previously chosen cell model
    # needed to adjust for the reshaping in the previous step
    cell$windows$NR_wmsh  <- with(cell, { windows$NR_wms * ploidy / windows$HMM })  
    cell$windows$NR_wmshf <- round(cell$windows$NR_wmsh[shapeFitI], 0) 
    cell$windows$NR_wmhfl <- excludeOutliers(cell$windows$NR_wmshf) # names intentionly "wrong" (no s), since used by modeling functions
    cell$windows$NR_wmhfld <- cell$windows$NR_wmhfl # use all windows, just solving one model
    cell <- setBestReplicationModel(cell, "reshaped")

    # adjust the states and state transition weights for the target replication model
    cell$states <- copy(compositeStates) # NAME	CN	NAR	N_ALLELES	CNC_P1	CNC_P2
    if(!cell$cellIsReplicating) {
        statesI <- cell$states[, NAR == 0]
        cell$states <- cell$states[statesI]
    }
    nStates <- nrow(cell$states)
    transProbs <- matrix(NA, nrow = nStates, ncol = nStates)
    for(fromI in 1:nStates) for(toI in 1:nStates){
        fromState <- cell$states[fromI]
        toState   <- cell$states[toI]
        transProbs[fromI, toI] <- if(fromI == toI) 0.99  # TODO: this should add up to 1 for each row/column
        else if(fromState$CN != toState$CN) 1e-8
        else 1e-1
    }

    # calculate the emission probabilities using all prior fitting as guidance
    theta <- cell$replicationModel$theta
    RPA <- predict(cell$gc_fit, gc_w)
    NR_wms <- round(cell$windows$NR_wms, 0)
    emissProbs <- sapply(cell$states$NAME, function(stateName){ # NAME	CN	NAR	N_ALLELES	CNC_P1	CNC_P2
        state <- cell$states[NAME == stateName]
        ER <- RPA * if(state$N_ALLELES == 0) 0.05 else state$N_ALLELES
        pNAR <- sapply(P_replicated_w, function(gci) P_replicated_wf[[gci]][[state$CN + 1]][state$NAR + 1])
        dnbinom(NR_wms, mu = ER, size = theta) * pNAR
    })

    # solve the HMM
    hmm <- new_hmmEPTable(log(emissProbs), transProbs = log(transProbs), keys = windows$chrom) 
    cell$windows$stateI <- keyedViterbi(hmm)
    cell$windows$NAR <- cell$states[cell$windows$stateI, NAR]
    plotWindows_counts(cell, "finalByRep", cell$windows$NR_wms, col = "useRepColor") 
    plotWindows_counts(cell, "finalByCn",  cell$windows$NR_wms, col = "useCnColor") 
    setHMMProfiles(cell)
}
#=====================================================================================

#=====================================================================================
# finish up by converting the composite HMM results to isolated CN and FAR profiles
#-------------------------------------------------------------------------------------
setHMMProfiles <- function(cell){
    # message("setHMMProfiles")
    windows <- windows[[cell$windowPower + 1]]
    cell$windows$HMM  <- cell$states[cell$windows$stateI, CN]
    cell$windows$NAR  <- cell$states[cell$windows$stateI, NAR]
    cell$windows$NA_  <- cell$states[cell$windows$stateI, N_ALLELES]
    cell$windows$NR_wmsh <- with(cell, { windows$NR_wms * ploidy / windows$HMM }) # correct NR toward ploidy, for analyzing replication
    cell$windows$NR_wmsr <- with(cell, { windows$NR_wms * windows$HMM / windows$NA_ }) # correct NR toward unreplicated, for analyzing CNVs
    cell$windows$CN      <- with(cell, { windows$NR_wmsr / windows$RPA }) # calculated copy number after accounting for replication
    cell$windows$FAR     <- with(cell, { (windows$NR_wmsh - windows$ER_ploidy) / windows$ER_ploidy }) # fraction replicated per window
    cell$fractionS <- with(cell$windows, { sum(NAR, na.rm =TRUE) / sum(CN, na.rm =TRUE) }) # FAR for the whole cell
    plotWindows_cn(cell, "isolated", windows[, gc_fraction], cell$windows$NR_wmsr, hmm = cell$windows$HMM)
    plotWindows_FAR(cell, "isolated")
    cell
}#=====================================================================================

#=====================================================================================
# main target function that launches the process!
#-------------------------------------------------------------------------------------
fitCell <- function(cell_id){ 
    # message("--------------------------------------------------")
    # message(cell_id)    
    plotNumber <<- 1 # so plots are ordered in lists as they were created
    cell <- list(
        cell_id = cell_id,
        modal_NA = env$PLOIDY, # subject to change if late S-phase
        ploidy   = env$PLOIDY  # never changed by this script, always env$PLOIDY
    )
    cell$minWindowCount <- getMinWindowCount(env$PLOIDY)
    cell$minWindowPower <- getMinWindowPower(cell_id, cell$minWindowCount)
    cell <- setCellWindows(cell) # cascades forward through all actions...
    for(key in c("CN_c", "NR_wmh", "NR_wmhf", "NR_wmhfl", "NR_wmhfld", "NR_wmhr", # remove data not needed for package
                "shape", "NR_wmshr", "NR_wmsh", "NR_wmshf", "NR_wmsr")) cell$windows[[key]] <- NULL
    windows <- as.list(rep(NA, env$MAX_WINDOW_POWER + 1))
    windows[[cell$windowPower + 1]] <- cell$windows
    cell$windows <- windows
    cell
}
#=====================================================================================


# List of 12
#  $ cell_id          : chr "3"
#  $ modal_NA         : num 4  # most frequent number of alleles per window
#  $ ploidy           : int 2  # how many copies of autosomes is typical
#  $ minWindowCount   : num 64 # least possible number of bins required to give good CNV calling/visualization
#  $ minWindowPower   : num 0  # log2(least number of bins to achieve minWindowCount)
#  $ windowPower      : num 2  # the optimized window power to give best CNV calling/visualization, filled in below
#  $ fractionS        : num 0.868 # extent of DNA replication, as calculated from the composite HMM
#  $ cellIsReplicating: logi TRUE # shortcut for interpreting replicationModel$modelType
#  $ windows          :List of 8 # indexed as windowPower + 1
#   ..$ : logi NA # provides empty windows sizes in case the app user wants to fill them in
#   ..$ : logi NA
#   ..$ :List of 24
#   .. ..$ NR_wm          : num [1:35614] 212 242 147 214 255 ... # the input counts with only mappability correction
#   .. ..$ ER_modal_NA    : num 284
#   .. ..$ ER_ploidy      : num 142
#   .. ..$ RPA            : num 70.9
#   .. ..$ normLagDiffQ   : Named num 0.171
#   .. .. ..- attr(*, "names")= chr "50%"
#   .. ..$ CN_c # DELETED
#   .. ..$ HMM            : int [1:35614] 2 2 2 2 2 2 2 2 2 2 ...
#   .. ..$ NR_wmh # DELETED
#   .. ..$ NR_wmhf # DELETED
#   .. ..$ NR_wmhfl # DELETED
#   .. ..$ NR_wmhfld # DELETED
#   .. ..$ NAR            : int [1:35614] 2 2 2 2 2 2 2 2 2 2 ...
#   .. ..$ NR_wmhr # DELETED
#   .. ..$ shape # DELETED
#   .. ..$ shapeCorrection: num [1:35614] 1.06 1.06 1.06 1.06 1.06 ...
#   .. ..$ NR_wms         : num [1:35614] 205 234 142 207 246 ...
#   .. ..$ NR_wmshr # DELETED
#   .. ..$ NR_wmsh # DELETED
#   .. ..$ NR_wmshf # DELETED
#   .. ..$ stateI         : num [1:35614] 6 6 6 6 6 6 6 6 6 6 ...
#   .. ..$ NA_            : int [1:35614] 4 4 4 4 4 4 4 4 4 4 ...
#   .. ..$ NR_wmsr # DELETED
#   .. ..$ CN             : num [1:35614] 1.44 1.65 1 1.46 1.74 ...
#   .. ..$ FAR            : num [1:35614] 0.443143 0.651339 0.000733 0.461978 0.73702 ...
#   ..$ : logi NA
#   ..$ : logi NA
#   ..$ : logi NA
#   ..$ : logi NA
#   ..$ : logi NA
#  $ gc_fit           :List of 8 # final GC bias model, derived from reshaped NR_wms vs. GC
#   ..$ nGcSteps     : num 100
#   ..$ gcIndexOffset: num 31
#   ..$ gcFractions  : num [1:35] 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 ...
#   ..$ minGcIndex   : num 32
#   ..$ maxGcIndex   : num 66
#   ..$ theta        : num [1:35] 38 23.6 23.7 27.8 28.6 ...
#   ..$ mu           : num [1:35] 58.1 61 63.7 66.1 68.4 ...
#   ..$ peak         : int [1:35] 56 58 61 63 65 68 69 71 67 67 ...
#   ..- attr(*, "class")= chr "nbinomCountsGC2"
#  $ replicationModel :List of 5 # final replication model, derived from reshaped NR_wms density
#   ..$ modelType      : chr "peakIsReplicated"
#   ..$ fractionS      : num 0.825 # the estimate from optim
#   ..$ ER_unreplicated: num 147
#   ..$ theta          : num 28.2
#   ..$ logLikelihood  : num -186921
#  $ states           :Classes ‘data.table’ and 'data.frame':     15 obs. of  6 variables: # cell-specific (e.g., (un)replicated)
#   ..$ NAME     : chr [1:15] "CN0_NAR0" "CN1_NAR0" "CN1_NAR1" "CN2_NAR0" ...
#   ..$ CN       : int [1:15] 0 1 1 2 2 2 3 3 3 3 ...
#   ..$ NAR      : int [1:15] 0 0 1 0 1 2 0 1 2 3 ...
#   ..$ N_ALLELES: int [1:15] 0 1 2 2 3 4 3 4 5 6 ...
#   ..$ CNC_P1   : int [1:15] -1 0 0 1 1 1 2 2 2 2 ...
#   ..$ CNC_P2   : int [1:15] -2 -1 -1 0 0 0 1 1 1 1 ...
#   ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "index")= int(0) 
#   .. ..- attr(*, "__NAME")= int(0) 

# # make a composite plot of cell model for QC purposes
# plotCellByWindow <- function(d, ylab, ylim = NULL){
#     plot(1:length(d), d, bty = "n",
#         pch = 19, cex = 0.4, 
#         col = rgb(0, 0, 0, 0.1),
#         xaxt = "n", xlab = NULL, ylab = ylab, ylim = ylim)
# }
# plotCellQC <- function(cell_id, cell){
#     png(
#         filename = paste(env$PLOT_PREFIX, cell_id, "qc", "png", sep = "."),
#         width  = 1.5 * 6, 
#         height = 1.35, 
#         units = "in", 
#         pointsize = 7,
#         bg = "white",
#         res = 96, # i.e., optimized for on-screen display
#         type = "cairo"
#     )
#     layout(matrix(c(c(1,1), rep(c(2,3), 5)), nrow = 2, ncol = 6))

#     # plot NR_map_w vs. gc_w
#     par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
#     plot(cell$gc_fit, cell$gc_w, cell$NR_map_w, cell$modal_NA, !cell$keep)

#     # plot NR_map_w vs. window index, i.e., pre-normalization
#     par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
#     plotCellByWindow(cell$NR_map_w, "# Reads", ylim = c(0, cell$ER_ploidy * 3))

#     # plot CN vs. window index, i.e., post-normalization
#     plotCellByWindow(cell$cn, "CN", ylim = c(0, 6))
#     abline(h = 0:4, col = "grey")
#     if(!is.null(cell$hmm)) lines(1:length(cell$hmm), cell$hmm, col = "red")

#     dev.off()
# }
