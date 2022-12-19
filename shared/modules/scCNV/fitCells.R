#=====================================================================================
# throughout this script and pipeline, the following naming conventions are used
#-------------------------------------------------------------------------------------
#   variable names:
#     NR  = Number of Reads, i.e., a read count (can be fractional due to fragment splitting across bins)
#     ER  = Expected Reads, i.e., the nominal read count that corresponds to a specific copy number and/or state
#     RPA = Reads Per Allele, i.e., the count expected for CN == 1, i.e., ER_ploidy / ploidy
#     CN  = Copy Number (also cn), PRIOR to any replication, can be fractional as NR / ER
#     HMM = Hidden Markov Model (also hmm), i.e., quantal CN values determined from a model of NR + ER
#     NA  = Number of Alleles, i.e., how many DNA copies are present considering BOTH CN and replication
#     NAR = Number of Alleles Replicated, thus, could be c(0,1,2) for CN==2, to yield NA = c(2,3,4)
#     FAR = Fraction of Alleles Replicated, i.e., mean(NAR) / ploidy, synonymous with fractionS
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
#     h = counts that have been rescaled to an HMM, e.g., NR * ploidy / hmm
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
#                 shape effects vary considerable with the amplification method
#   - "fractionS" is the total fraction of alleles across the genome that have been replicated (FAR), e.g.,
#                 at fractionS==0.5, half of the genome's bases (not windows!) have been replicated
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
saveCellPlot <- function(cell, name, fn, width = 2, height = 2){ # save a plot for future assembly and app
    # TODO: update plotDir to proper output directory
    plotDir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV/plots"
    plotDir <- file.path(plotDir, cell$cell_id)
    if(!dir.exists(plotDir)) dir.create(plotDir)
    filename <- paste(cell$cell_id, letters[plotNumber], name, cell$windowPower, "png", sep = ".")
    pngFile <- file.path(plotDir, filename)
    png(pngFile, width = width, height = height, units = "in", pointsize = 7, res = 300, type = "cairo")
    fn()
    dev.off()
    plotNumber <<- plotNumber + 1
}
getRepColor <- function(cell, default = NULL){ # after replication fitting, color points by replication state
    cols <- c(rgb(0, 0, 1, 0.1), rgb(0, 1, 0, 0.1), rgb(1, 0, 0, 0.1))  
    if(is.null(cell$replicationModel)) {
        if(is.null(default)) defaultPointColor else default
    } else cols[cell$replicationModel$NAR + 1]
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
        plot(NA, NA, xlim = coord$xlim, ylim = c(0, cell$ER_modal_CN * 3),
             xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "# of Reads")
        addChromLabels(coord)
        if(is.null(cell$RPA)) abline(h = cell$RPA * 0:maxModelCN, col = "grey")        
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
        abline(h = 0:maxModelCN, col = "grey")
        points(coord$x, cn, pch = 19, cex = 0.3, col = col)
        if(!is.null(hmm)) lines(coord$x, hmm, col = "red3")
    }, width = 6)
}
plotCellFit <- function(cell, model, bestModelByType){ # distributions of read counts with best-fitted replication models
    saveCellPlot(cell, "repFit", function(){
        NR_a <- round(cell$NR_wmhfl, 0) # here, a is for "actual"
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
        abline(v = cell$ER_modal_CN * c(1/2, 1, 2))
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
# places ~96% of windows within (modal_CN +/- 0.5) * replicationFactor==[1,2]
# sapply(1:4, getMinWindowCount) => 16 64 144 256; thus, usually 64 reads/window for interphase diploid cells
getMinWindowCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
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
getCellWindows <- function(cell, windowPower){ # parse a cell's windows at a given windowPower
    windowSize <- 2 ** windowPower
    windows <- windows[[windowPower + 1]]
    mappability <- windows[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_wm <- unname(unlist(sapply(constants$chrom, function(chrom){
        collapseVector(
            raw_counts[[cell$cell_id]][rowRanges$chrom == chrom], 
            windowSize
        ) / mappability[windows$chrom == chrom]
    })))
    NR_wmf <- excludeOutliers(NR_wm, min = 1) # thus, chrom-to-chrom, baseline-to-CNV, and rare bad windows won't have undue influence
    ER_modal_CN <- peakValue(NR_wmf) # ER_modal_CN may NOT be the same as ER_ploidy 
    list( # one of these is chosen to become the working cell object
        cell_id = cell$cell_id,
        windowPower = windowPower,
        NR_wm       = NR_wm,        
        modal_CN    = env$PLOIDY, # subject to change if late S-phase
        ploidy      = env$PLOIDY,
        ER_modal_CN = ER_modal_CN,
        ER_ploidy   = ER_modal_CN,
        RPA         = ER_modal_CN / env$PLOIDY, # later we may learn RPA is 2-fold off for replicating cells,
        normLagDiffQ = getNormLagDiffQ(NR_wmf)
    )
}
setCellWindows <- function(cell){ # establish the optimal window power for a cell
    message("setCellWindows")
    x <- list() 
    windowPower <- cell$minWindowPower
    x[[windowPower + 1]] <- getCellWindows(cell, windowPower)
    while(x[[length(x)]]$normLagDiffQ > normLagDiffQ_threshold && # stop when we get a nice, appropriately tight distribution
          windowPower < env$MAX_WINDOW_POWER){
        windowPower <- windowPower + 1
        x[[windowPower + 1]] <- getCellWindows(cell, windowPower)
    }
    x
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
    cell$hmm <- viterbi(
        cell$gc_fit, 
        NR_wmf, # note, could simply be NR_wm if windowI == TRUE
        gc_wf, 
        maxCN = maxModelCN, 
        chroms = windows[windowI == TRUE, chrom], 
        asRle = FALSE
    )$cn
    plotWindows_cn(cell, name, gc_wf, NR_wmf, windowI, cell$hmm)
    cell
}
#=====================================================================================

#=====================================================================================
# use a replication-unaware algorithm that squashes the replication GC bias effect
# to obtain an initial copy number estimates across the genome
#-------------------------------------------------------------------------------------
maxModelCN <- 5
CNs <- 0:maxModelCN # CN==5 is not trustworthy, the true value could be (much) higher
setChromCN <- function(cell){ # cell$windows called just cell here for readability, here and below
    message("setChromCN")
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # determine a first crude aneuploidy copy number estimate for each whole chromosome
    NR_wmf <- excludeOutliers(cell$NR_wm, min = 2)
    if(length(NR_wmf) > 1e4) NR_wmf <- sample(NR_wmf, 1e4)
    RPA_wmf <- NR_wmf / cell$ploidy
    NR_wmf_cn <- lapply(CNs, "*", RPA_wmf)
    cell$CN_c <- lapply(chroms, function(chrom){
        NR_wmcl <- excludeOutliers(cell$NR_wm[windows$chrom == chrom], min = 2)
        if(length(NR_wmcl) <= 1 || median(NR_wmcl) <= 1) return(0)
        p <- sapply(NR_wmf_cn, function(NR_g) wilcox.test(NR_g, NR_wmcl, exact = FALSE)$p.value)
        CNs[which.max(p)] # the copy number that gives the best distribution match between chrom and genome
    })
    names(cell$CN_c) <- chroms
    # cell$chromCnIsTrustworthy <- cell$CN_c[windows$chrom] %between% c(1, maxModelCN - 1)

    # use the crude estimates to select chromosomes for performing a proper GC bias fit and CN HMM
    chromFitI <- unlist(cell$CN_c) == cell$ploidy
    col_w <- sapply(windows$chrom, function(chrom) {
        if(cell$CN_c[[chrom]] == cell$ploidy) defaultPointColor else rgb(1, 0, 0, 0.1)
    })
    cell <- fitGcBias(cell, "input", cell$NR_wm, chroms = chroms[chromFitI], col_w = col_w)
    cell <- solveSimpleHMM(cell, "squashed", cell$NR_wm)
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
getModelER_unreplicated <- function(modelType, cell){ # set read count expectations (ER) for a combination of modelType and cell$ER_modal_CN
    switch(
        modelType,
        notReplicating     = cell$ER_modal_CN,
        peakIsUnreplicated = cell$ER_modal_CN,    # the modal window count corresponds to unreplicated DNA
        peakIsReplicated   = cell$ER_modal_CN / 2 # the modal window count corresponds to   replicated DNA
    )
}
fillModelLikelihoods <- function(modelType, fractionS, cell){

    # use a rich subset of the data to fit the candidate model to the data as best as possible
    # models are defined by replication peak type and fractionS increments
    ER_unreplicated <- getModelER_unreplicated(modelType, cell)
    optFn <- function(par){ 
        LL <- -getModelLogLikelihood(
            dnbinom(cell$NR_wmhfld, mu = par[1],       size = par[2]), 
            dnbinom(cell$NR_wmhfld, mu = par[1] * 1.5, size = par[2]),
            dnbinom(cell$NR_wmhfld, mu = par[1] * 2,   size = par[2]),
            fractionS
        )
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
    par[3] <- getModelLogLikelihood(
        dnbinom(cell$NR_wmhfl, mu = par[1],       size = par[2]), 
        dnbinom(cell$NR_wmhfl, mu = par[1] * 1.5, size = par[2]), 
        dnbinom(cell$NR_wmhfl, mu = par[1] * 2,   size = par[2]),
        fractionS
    )
    par
}
getRepEmissProbs <- function(cell){
    model <- cell$replicationModel
    weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, model$fractionS)
    log(as.matrix(data.table( # here, model$ER_unreplicated has been fitted to reflect mu, not peak/mode
        dnbinom(cell$NR_wmhf, mu = model$ER_unreplicated,       size = model$theta) * weights[1], 
        dnbinom(cell$NR_wmhf, mu = model$ER_unreplicated * 1.5, size = model$theta) * weights[2], 
        dnbinom(cell$NR_wmhf, mu = model$ER_unreplicated * 2,   size = model$theta) * weights[3]
    )))
}
#=====================================================================================

#=====================================================================================
# fit a cell's data to establish the initial parameters of it's replication model
#-------------------------------------------------------------------------------------
setReplicationModel <- function(cell){
    message("setReplicationModel")
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # use CNV-corrected counts to establish cell-specific likelihoods for replication models
    cell$NR_wmh <- with(cell, { NR_wm * ploidy / hmm }) 
    repFitI <- with(cell, { hmm > 0 & hmm < maxModelCN & !is.na(NR_wmh) }) 
    cell$NR_wmhf <- round(cell$NR_wmh[repFitI], 0)
    cell$NR_wmhfl <- excludeOutliers(cell$NR_wmhf)
    cell$NR_wmhfld <- with(cell, { if(length(NR_wmhfl) > 1e4) sample(NR_wmhfl, 1e4) else NR_wmhfl }) 

    # optimize theta for three distinct types of models, at multiple fixed fractionS values
    # and ER_(un)replicated values determine by model type
    models <- copy(modelValues)
    models <- cbind(models, t(models[, mapply(function(modelType, fractionS){
        fillModelLikelihoods(modelType, fractionS, cell)
    }, modelType, fractionS)]))
    setnames(models, c("modelType", "fractionS", "ER_unreplicated", "theta", "logLikelihood"))

    # pick the best model type, with its best parameters
    cell$replicationModel <- as.list(models[which.max(logLikelihood)])
    bestModelByType <- models[, {
        i <- which.max(logLikelihood)
        .SD[i]
    }, by = "modelType"]
    plotCellFit(cell, cell$replicationModel, bestModelByType) # based on NR_wmhfl

    # calculate derived model parameters
    peakIsReplicated <- cell$replicationModel$modelType == "peakIsReplicated"
    cell$modal_CN  <- cell$ploidy * (if(peakIsReplicated) 2 else 1)
    cell$ER_ploidy <- cell$ER_modal_CN / (if(peakIsReplicated) 2 else 1)
    cell$RPA <- cell$ER_modal_CN / cell$modal_CN # all not yet subjected to gc fit
    cell$cellIsReplicating <- cell$replicationModel$modelType != "notReplicating"

    # solve an HMM to establish the (un)replicated genome spans for adjusting window NR prior to GC fit
    message("solve rep HMM")
    emissProbs <- getRepEmissProbs(cell)
    hmm <- new_hmmEPTable(emissProbs, transProb = 1e-1, keys = windows$chrom[repFitI])
    cell$replicationModel$NAR <- {
        NAR <- keyedViterbi(hmm) - 1 # based on NR_wmhf
        round(NAR * cell$hmm[repFitI] / cell$ploidy, 0) # rescale NAR to hmm
    }
    plotWindows_counts(cell, "repHMM", cell$NR_wm[repFitI], windowI = repFitI, col = "useRepColor") 
    setChromShapes(cell, repFitI) # pass the baton
}
#=====================================================================================

#=====================================================================================
# apply a further correction to replication squashed data that accounts for the shape
# of individual chromosome coverage profiles, which results from amplification artifacts
#-------------------------------------------------------------------------------------
setChromShapes <- function(cell, shapeFitI){ # shapeFitI is the same as repFitI
    message("setChromShapes")
    windows <- windows[[cell$windowPower + 1]]

    # fit a quadratic to each chromosome and use to normalize counts, i.e., stabilize the NR baseline
    NAR <- rep(NA, nrow(windows)) # begin with a list of all windows, only some of which have NAR data
    NAR[shapeFitI] <- cell$replicationModel$NAR
    cell$NR_wmhr <- with(cell, { NR_wm * modal_CN / (hmm + NAR) }) # thus, correct counts toward ploidy
    cell$shape <- windows[, { # shape is relative to ER_modal_CN
        sfI  <- shapeFitI[.I] # model is created from filtered bins
        if(sum(sfI) < 10) rep(1, .N) else { # insufficient points to fit
            i_f  <- i[sfI]
            i2_f <- i2[sfI]
            NR_wmhrf <- cell$NR_wmhr[.I][sfI]
            if(all(is.na(NR_wmhrf))) rep(1, .N)
            fit <- lm(NR_wmhrf ~ i_f + i2_f)
            predict(fit, newdata = data.frame(i_f = i, i2_f = i2)) # final shape prediction occurs on ALL bins 
        }
    }, by = "chrom"][[2]]
    cell$shapeCorrection <- cell$shape / cell$ER_modal_CN        
    cell$NR_wms <- with(cell, { NR_wm / shapeCorrection })
    cell$NR_wms <- with(cell, { NR_wms * sum(NR_wm, na.rm = TRUE) / sum(NR_wms, na.rm = TRUE) }) # maintain total cell weights
    plotWindows_counts(cell, "shapes", cell$NR_wmhr[shapeFitI],
                       windowI = shapeFitI, shape = cell$shape[shapeFitI])

    # re-solve the GC bias fit using the reshaped NR values
    cell$NR_wmshr <- with(cell, { NR_wms * ploidy / (hmm + NAR) }) # thus, correct counts toward ploidy
    col_w <- ifelse(shapeFitI, defaultPointColor, rgb(1, 0, 0, 0.1))
    cell <- fitGcBias(cell, "reshaped", cell$NR_wmshr, windowI = shapeFitI, col_w = col_w)
    plotWindows_counts(cell, "reshaped", cell$NR_wms[shapeFitI], windowI = shapeFitI, col = "useRepColor")  
    solveCompositeHmm(cell, shapeFitI) # pass the baton
}
#=====================================================================================

#=====================================================================================
# support functions for replication state profiles as a function of fractionS + window GC
# i.e., that add a (P_unreplicated | GC) weight to the composite HMM for replicating cells
# derived from intial fit of each individual cell
#-------------------------------------------------------------------------------------
nGcSteps <- 50
gcIndices  <- 1:nGcSteps
P_rep_fs_gc <- list() # for collecting replicating cell profiles for later aggregation
getP_unreplicated_wf <- function(cell, gc_wf) with(cell, {
    gci_wf <- round(gc_wf * nGcSteps, 0) 
    pmax(0.01, sapply(gcIndices, function(gci){ # all bins have at least somelikelihood of being called replicating
        NAR <- replicationModel$NAR[gci_wf == gci]
        if(length(NAR) == 0) return(NA)
        agg <- aggregate(NAR, list(NAR), length)
        n0 <- agg[agg[[1]] == 0, 2]
        n0 <- if(length(n0) == 0) 0 else n0
        n0 / sum(agg[[2]]) # thus, one P_unreplicated value per fraction GC bin
    }))
})
getP_unreplicated_w <- function(cell, gc_w)  with(cell$replicationModel, {
    gci_w <- round(gc_w * nGcSteps, 0)
    is <- 1:length(gci_w)    
    nonNaGci <- which(!is.na(P_unreplicated_wf))  
    sapply(is, function(i){ # thus, one P_unreplicated per window
        gci <- gci_w[i]
        p <- P_unreplicated_wf[gci]
        if(length(p) > 0 && !is.na(p)) return(p)
        delta <- abs(nonNaGci - gci)
        gci <- nonNaGci[which(delta == min(delta))[1]]
        P_unreplicated_wf[gci]
    })
})
commitReplicationProfile <- function(cell, gc_wf) with(cell, {
    if(!cellIsReplicating) return(NULL)
    gcIndex <- round(gc_wf * nGcSteps, 0) 
    P_rep_fs_gc[[cell_id]] <<- list(
        cell_id = cell_id,
        replicationModel = replicationModel,
        profile = sapply(gcIndices, function(gci){
            nAllelesReplicated <- replicationModel$nAllelesReplicated[gcIndex == gci]
            if(length(nAllelesReplicated) == 0) return(c(0, NA, NA, NA))
            agg <- aggregate(nAllelesReplicated, list(nAllelesReplicated), length)
            N <- sum(agg[[2]])
            c(N, sapply(0:2, function(nar) {
                nnar <- agg[agg[[1]] == nar, 2]
                if(length(nnar) == 0) 0 else nnar
            }) / N)
        })
    )  
})
plotReplicationProfiles <- function(){
    values <- c("N", "P_NAR0", "P_NAR1", "P_NAR2")
    fractionS_int <- as.integer(seq(0.05, 0.95, 0.025) * 1000)
    colfunc <- colorRampPalette(c("blue", "red"))
    colors <- colfunc(length(fractionS_int))
    weights <- sapply(P_rep_fs_gc, function(cell) cell$profile[1,])
    for(j in 1:length(values)){ # one plot per value type
        message(values[j])
        saveDevPlot(paste("P_rep_fs_gc", values[j], sep = "."), function(){
            x  <- gcIndices
            ys <- sapply(P_rep_fs_gc, function(cell) cell$profile[j,])
            plot(
                NA, NA, 
                xlim = c(0.3, 0.7), ylim = range(ys, na.rm = TRUE), 
                xlab = "Fraction GC", ylab = values[j]
            ) 
            for(k in 1:ncol(ys)) { # once trance per cell
                fractionS <- 
                P_rep_fs_gc[[k]]$replicationModel$fractionS
                col <- colors[fractionS_int == as.integer(fractionS * 1000)]
                y <- ys[, k]
                lines(x / nGcSteps, y, col = col)  
                points(x / nGcSteps, y, col = col)               
            }
        })
    }
}
#=====================================================================================

#=====================================================================================
# use the replication model and reshaped GC bias fit to solve a multi-state HMM with CN and replication components
#-------------------------------------------------------------------------------------
solveCompositeHmm <- function(cell, shapeFitI){
    windows <- windows[[cell$windowPower + 1]]
    gc_w <- windows[, gc_fraction]
    gc_wf <- gc_w[shapeFitI]
    cell$replicationModel$P_unreplicated_wf <- getP_unreplicated_wf(cell, gc_wf)
    P_unreplicated_w <- getP_unreplicated_w(cell, gc_w)

    str(cell$replicationModel)
    print(cell$replicationModel$P_unreplicated_wf)
    str(P_unreplicated_w)

    stop("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

    # thus, use 
    #   cell$NR_wms to reflect chromosome shape, i.e., amplification skew
    #   cell$gc_fit to determine ER per state as gc_fit(gc) * (cn + nar = na)
    #   P_unreplicated_wf to weight windows for replication
    # NONE of these are deterministic for a window! rather, they are trends determined from aggregated data
    # thus, the model output is not dictated or forced at any window by prior preparative work

    p_state <- dnbinom(cell$NR_wms, mu = ER, size = theta) * P_unreplicated_w

    # question: re-calculate theta after reshaping? likely yes!
    # only need to solve for the actual modelType, with fractionS similar to established value
    # main goal is to establish theta based on NR_wms

    # load transition weight table, may need to optimize this
    # the key to proper output will be preventing NAR to NAR transitions as being CNVs
    # thus, in mid-S, NAR transitions should be more likely than CNV transitions
}
#=====================================================================================

# fit a GC bias correction to the normalized window data and solve for CNVs by HMM
fitCellGcHmm <- function(cell){
    windows <- windows[[cell$windowPower + 1]]
    repFitI <- !is.na(cell$replicationModel$nAllelesReplicated)
    NR_map_wcr <- NR_map_wc / cell$replicationModel$nAllelesReplicated
    
    with(cell, { NR_map_wc / ifelse(replicationModel$nAlleleReplicated, 2, 1) }) # thus, these are adjusted toward ploidy
    NR_map_wcrf <- NR_map_wcr[repFitI]    

    # use aneuploidy + replication-corrected counts to fit a cell's GC bias using the negative binomial distribution
    gc_wf <- windows[repFitI == TRUE, gc_fraction]
    cell$fit <- new_nbinomCountsGC2(NR_map_wcrf, gc_wf, binCN = cell$ploidy)



    # solve a final CN estimate for all windows
    cell$gc_w <- windows[, gc_fraction] 
    cell$ER_gc <- predict(cell$fit, cell$gc_w, type = 'adjustedPeak') * cell$ploidy 
    cell$cn <- cell$NR_map_w / cell$ER_gc * cell$ploidy # thus, in late S modal_CN windows will have cn == 4
    cell$theta <- predict(cell$fit, cell$gc_w, type = 'theta')

    # ALTERNATIVE: solve a multi-state model of CN +/- replication, instead of correcting back and forth
    
    # correct for replication prior to running CNV HMM
    cell$NR_map_wr <- cell$NR_map_w / ifelse(cell$replicationModel$windowIsReplicated, 2, 1)
    cell$hmm <- viterbi(cell$fit, cell$NR_map_wr, cell$gc_w, asRle = FALSE, 
                        chroms = windows$chrom, transProb = env$TRANSITION_PROBABILITY)$cn
    cell$percentile <- cumprob(cell$fit, cell$NR_map_wr, cell$gc_w, binCN = cell$hmm)    

    # correct for CNVs prior to re-running replication HMM
    # cell$NR_map_wc <- cell$NR_map_w * cell$ploidy / cell$hmm 
    # cell$NR_map_wc[cell$hmm == 0] <- NA
    # TODO: this is returning cn, not logical replication
    # cell$windowIsReplicated <- if(!cell$replicationModel$cellIsReplicating) rep(cell$ploidy, length(cell$NR_map_wc))
    #                            else viterbi(cell$fit, cell$NR_map_wc, cell$gc_w, asRle = FALSE, 
    #                                          chroms = windows$chrom, transProb = env$TRANSITION_PROBABILITY)$cn

    # repeat fit and HMMs once more to refine?
    cell
}

fitCell <- function(cell_id, stage = "extract"){ 
    plotNumber <<- 1
    cell <- list(cell_id = cell_id)

    message("--------------------------------------------------")
    message(cell_id)

    cell$minWindowCount <- getMinWindowCount(env$PLOIDY)
    cell$minWindowPower <- getMinWindowPower(cell_id, cell$minWindowCount)
    cell$windows <- setCellWindows(cell) # will contain some subset of windowPowers
    cell$windowPower <- length(cell$windows) - 1
    windowPowerIndex <- cell$windowPower + 1    

    print(cell$windowPower)  


    cell$windows[[windowPowerIndex]] <- setChromCN(cell$windows[[windowPowerIndex]]) # cascades forward to next actions...
    #cell$windows[[windowPowerIndex]] <- setReplicationModel(cell$windows[[windowPowerIndex]])
    # cell$windows[[windowPowerIndex]] <- fitCellGcHmm(cell$windows[[windowPowerIndex]])    


    # str(minWindowCount)
    # str(minWindowPower)
    # str(cell)

    return(cell)

    stop("xxxxxxx")

    # TODO: run this twice
    # first time, use crude CN correction
    # second time, pass the CN obtained from HMM in first round
    # also pass the GC fit from 1st round to the 2nd round?
    # consider also resetting window power using GC fit 2nd time?

    # espeically if windowPower setting is unreliable, consider doing this for several windowPowers per cell? min:default

    cell <- fitCellGcHmm(cell)    

    ci = cell_id
    print(colData[cell_id == ci])
    # print(minWindowCount)
    # print(minWindowPower)
    str(cell)

    # TMP_DIR <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV"
    # nBins <- length(cell$NR_map_w)
    # col <- if(cell$replicationModel$cellIsReplicating)

    # pngFile <- file.path(TMP_DIR, paste0("NR_map_w", ".png"))
    # png(pngFile, width = 6, height = 3, units = "in", pointsize = 8, res = 300, type = "cairo")
    # plot(1:nBins, cell$NR_map_w, pch = 16, cex = 0.25, col = col)
    # dev.off()

    stop("XXXXXXXXXXXXXXXXXXXXXXXXXX")

}

# make a composite plot of cell model for QC purposes
plotCellByWindow <- function(d, ylab, ylim = NULL){
    plot(1:length(d), d, bty = "n",
        pch = 19, cex = 0.4, 
        col = rgb(0, 0, 0, 0.1),
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

    # plot NR_map_w vs. gc_w
    par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
    plot(cell$gc_fit, cell$gc_w, cell$NR_map_w, cell$modal_CN, !cell$keep)

    # plot NR_map_w vs. window index, i.e., pre-normalization
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plotCellByWindow(cell$NR_map_w, "# Reads", ylim = c(0, cell$ER_ploidy * 3))

    # plot CN vs. window index, i.e., post-normalization
    plotCellByWindow(cell$cn, "CN", ylim = c(0, 6))
    abline(h = 0:4, col = "grey")
    if(!is.null(cell$hmm)) lines(1:length(cell$hmm), cell$hmm, col = "red")

    dev.off()
}
