# determine the least number of reads per window for a cell's data to support a robust HMM
# places ~96% of windows within (modal_CN +/- 0.5) * replicationFactor==[1,2]
# sapply(1:4, getMinWindowCount) => 16 64 144 256; thus, usually 64 reads/window for interphase diploid cells
getMinWindowCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
    (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# further require that the majority of windows have non-zero count, with allowance for rare homozygous losses and chrY
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
    NR_raw_b <- raw_counts[[cell_id]][bins]
    NR_avg_b <- mean(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_avg_b) || NR_avg_b == 0) return(env$MAX_WINDOW_POWER)
    windowSize <- ceiling(minWindowCount / NR_avg_b) # in number of bins, not bp
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
peakValue <- function(x){
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
normLagDiffQ_quantile  <- 0.5 # 0.95 # determined empirically
normLagDiffQ_threshold <- 0.175 # 0.5  # determined empirically
getNormLagDiffQ <- function(NR){ # 95%ile of the magnitude of the window-to-window count difference, normalized to local mean count
    # TODO: exclude outliers?
    lagDiff <- abs(diff(NR))
    lagMean <- (head(NR, -1) + tail(NR, -1)) / 2
    quantile(lagDiff / lagMean, normLagDiffQ_quantile, na.rm = TRUE)
}
getCellWindows <- function(cell, windowPower){ # parse a cell's windows at a given windowPower
    windowSize <- 2 ** windowPower
    windows <- windows[[windowPower + 1]]
    mappability <- windows[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_map_w <- unname(unlist(sapply(constants$chrom, function(chrom){
        collapseVector(
            raw_counts[[cell$cell_id]][rowRanges$chrom == chrom], 
            windowSize
        ) / mappability[windows$chrom == chrom]
    })))
    NR_map_wf <- excludeOutliers(NR_map_w, min = 1) # thus, chrom-to-chrom, baseline-to-CNV, and rare bad windows won't have undue influence
    ER_modal_CN <- peakValue(NR_map_wf) # ER_modal_CN may NOT be the same as ER_ploidy
    readsPerAllele <- ER_modal_CN / env$PLOIDY # later we may learn readsPerAllele is 2-fold off for replicating cells
    list( # one of these is chosen to become the working cell object
        cell_id = cell$cell_id,
        windowPower = windowPower,
        NR_map_w = NR_map_w,        
        modal_CN = env$PLOIDY, # subject to change if late S-phase
        ploidy   = env$PLOIDY,
        ER_modal_CN = ER_modal_CN,
        ER_ploidy   = ER_modal_CN,
        readsPerAllele = readsPerAllele,
        normLagDiffQ = getNormLagDiffQ(NR_map_wf)
    )
}
setCellWindows <- function(cell){
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

# determine a first crude aneuploidy copy number estimate for each whole chromosome
# small CNVs will persist but don't change the fitting outcome much
maxModelCN <- 5
cns <- 0:maxModelCN # CN==5 is not trustworthy, the true value could be (much) higher
setChromCN <- function(cell){ # cell$windows called just cell here for readability, here and below
    message("setChromCN")
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]
    NR_map_wf <- excludeOutliers(cell$NR_map_w, min = 2)
    if(length(NR_map_wf) > 1e4) NR_map_wf <- sample(NR_map_wf, 1e4)
    RPA_wf <- NR_map_wf / cell$ploidy
    NR_wf_cn <- lapply(cns, "*", RPA_wf)
    cell$chromCn <- lapply(chroms, function(chrom){
        NR_chrom <- excludeOutliers(cell$NR_map_w[windows$chrom == chrom], min = 2)
        if(length(NR_chrom) <= 1 || median(NR_chrom) <= 1) return(0)
        p <- sapply(NR_wf_cn, function(NR_g) wilcox.test(NR_g, NR_chrom, exact = FALSE)$p.value)
        cns[which.max(p)] # the copy number that gives the best distribution match between chrom and genome
    })
    names(cell$chromCn) <- chroms
    cell$chromCnIsTrustworthy <- cell$chromCn[windows$chrom] %between% c(1, maxModelCN - 1)
    setReplicationModel(cell) # pass the baton to the next encapsulated action
}

# establish the common parameters of a set of replication models to be applied to individual cells
modelTypes <- list(
    notReplicating     = 0,                      # allowable fractionsS values by peak type
    peakIsUnreplicated = seq(0.05, 0.75, 0.025), # where fractionS is the total fractional content of replicated alleles
    peakIsReplicated   = seq(0.25, 0.95, 0.025)
)
modelValues <- do.call(rbind, lapply(names(modelTypes), function(modelType){
    data.table(  # pre-assemble a table of all models
        modelType = modelType,
        fractionS = modelTypes[[modelType]]
    )
}))
modelTypeColors <- list( # for debugging plots
    notReplicating      = "blue",
    peakIsUnreplicated  = "green",
    peakIsReplicated    = "red"
)
fractionSFile <- file.path(env$ACTION_DIR, "fractionS_table.rds")
fractionSLookup <- readRDS(fractionSFile) # for converting overall fractionS to fraction of 0, 1, and 2 alleles replicated

# support functions for fitting replication models
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
            dnbinom(cell$NR_map_wcfld, mu = par[1],       size = par[2]), 
            dnbinom(cell$NR_map_wcfld, mu = par[1] * 1.5, size = par[2]),
            dnbinom(cell$NR_map_wcfld, mu = par[1] * 2,   size = par[2]),
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
        dnbinom(cell$NR_map_wcfl, mu = par[1],       size = par[2]), 
        dnbinom(cell$NR_map_wcfl, mu = par[1] * 1.5, size = par[2]), 
        dnbinom(cell$NR_map_wcfl, mu = par[1] * 2,   size = par[2]),
        fractionS
    )
    par
}
getRepEmissProbs <- function(cell){
    model <- cell$replicationModel
    weights <- fractionSLookup$getRepStateProbs(fractionSLookup, env$PLOIDY, model$fractionS)
    log(as.matrix(data.table( # here, model$ER_unreplicated has been fitted to reflect mu, not peak/mode
        dnbinom(cell$NR_map_wcf, mu = model$ER_unreplicated,       size = model$theta) * weights[1], 
        dnbinom(cell$NR_map_wcf, mu = model$ER_unreplicated * 1.5, size = model$theta) * weights[2], 
        dnbinom(cell$NR_map_wcf, mu = model$ER_unreplicated * 2,   size = model$theta) * weights[3]
    )))
}

# fit a cell's data to establish the parameters of it's replication model
setReplicationModel <- function(cell){
    message("setReplicationModel")
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # use aneuploidy-corrected counts to establish cell-specific likelihoods for replication models
    cell$NR_map_wc <- with(cell, { unname(unlist(sapply(chroms, function(chrom){ # rescale aneuploid chromosomes (including chrX) to ploidy
        NR_map_w[windows$chrom == chrom] * ploidy / chromCn[[chrom]]
    }))) })
    repFitI <- with(cell, { chromCnIsTrustworthy & !is.na(NR_map_wc) })
    cell$NR_map_wcf <- round(cell$NR_map_wc[repFitI], 0)
    cell$NR_map_wcfl <- excludeOutliers(cell$NR_map_wcf)
    cell$NR_map_wcfld <- if(length(cell$NR_map_wcfl) > 1e4) sample(cell$NR_map_wcfl, 1e4) else cell$NR_map_wcfl

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
    plotCellFit(cell, cell$replicationModel, bestModelByType)

    # calculate derived model parameters
    peakIsReplicated <- cell$replicationModel$modelType == "peakIsReplicated"
    cell$modal_CN  <- cell$ploidy * (if(peakIsReplicated) 2 else 1)
    cell$ER_ploidy <- cell$ER_modal_CN / (if(peakIsReplicated) 2 else 1)
    cell$readsPerAllele <- cell$ER_modal_CN / cell$modal_CN # all not yet subjected to gc fit

    # solve an HMM to establish the (un)replicated genome spans for adjusting window NR prior to GC git
    message("solve rep HMM")
    emissProbs <- getRepEmissProbs(cell)
    hmm <- new_hmmEPTable(emissProbs, transProb = 1e-1, keys = windows$chrom[repFitI])
    cell$replicationModel$nAllelesReplicated <- keyedViterbi(hmm) - 1

    # use aneuploidy + replication-corrected counts to fit a cell's GC bias using the negative binomial distribution
    repStateCorrection <- c(1, 1.5, 2)    
    repStateCn <- repStateCorrection * cell$ploidy
    gc_wf <- windows[repFitI == TRUE, gc_fraction]
    NR_map_wcfr <- with(cell, { # thus, window counts are adjusted toward ploidy
        NR_map_wcf / repStateCorrection[replicationModel$nAllelesReplicated + 1]
    })

    # probably need to fit gc_w (not f) setting NR to NA when filtered, to ensure a complete model?
    cell$gc_fit <- new_nbinomCountsGC2(NR_map_wcfr, gc_wf, binCN = cell$ploidy) # since corrected to ploidy

    plotGcBias("gc_uncorrected", cell, gc_wf, cell$NR_map_wcf)
    plotGcBias("gc_corrected",   cell, gc_wf, NR_map_wcfr)

    plotGenomeWindows(cell, cell$NR_map_wcf)
    plotGenomeWindows_cn(cell, gc_wf, cell$NR_map_wcf)

    return(cell)

    fitCellGcHmm(cell) # pass the baton to the next encapsulated action
}

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

saveDevPlot <- function(name, fn, width = 6){
    TMP_DIR <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV/plots"
    pngFile <- file.path(TMP_DIR, paste0(name, ".png"))
    png(pngFile, width = width, height = 3, units = "in", pointsize = 8, res = 300, type = "cairo")
    fn()
    dev.off()    
}
getAgg <- function(NR){
    x <- round(NR, 0)
    agg <- aggregate(x, list(x), length)
    names(agg) <- c("NR", "N")
    agg$freq <- agg$N / sum(agg$N)
    agg
}
plotCellFit <- function(cell, model, bestModelByType){
    saveDevPlot(paste(cell$cell_id, "plotCellFit", cell$windowPower, sep = "."), function(){
        message("plotCellFit")
        NR_a <- round(cell$NR_map_wcfl, 0)
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
    })
}
plotGcBias <- function(name, cell, gc_wf, NR_map_wcfr){
    saveDevPlot(paste(cell$cell_id, name, cell$windowPower, sep = "."), function(){
        message(name)
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        cols <- c(rgb(0, 0, 1, 0.1), rgb(0, 1, 0, 0.1), rgb(1, 0, 0, 0.1))    
        col <- cols[cell$replicationModel$nAllelesReplicated + 1]
        plot(cell$gc_fit, gc_wf, NR_map_wcfr, cell$ploidy, col = col)
    }, width = 3)
}
plotGenomeWindows <- function(cell, NR_map_wcf){
    saveDevPlot(paste(cell$cell_id, "genome_windows", cell$windowPower, sep = "."), function(){
        message("plotGenomeWindows")
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        cols <- c(rgb(0, 0, 1, 0.1), rgb(0, 1, 0, 0.1), rgb(1, 0, 0, 0.1))    
        col <- cols[cell$replicationModel$nAllelesReplicated + 1]
        plot(1:length(NR_map_wcf), NR_map_wcf, pch = 19, cex = 0.4, col = col)
        abline(h = cell$readsPerAllele * 0:6)
    })
}
plotGenomeWindows_cn <- function(cell, gc_wf, NR_map_wcf){
    saveDevPlot(paste(cell$cell_id, "plotGenomeWindows_cn", cell$windowPower, sep = "."), function(){
        message("plotGenomeWindows_cn")
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        cols <- c(rgb(0, 0, 1, 0.1), rgb(0, 1, 0, 0.1), rgb(1, 0, 0, 0.1))    
        col <- cols[cell$replicationModel$nAllelesReplicated + 1]

        ER_gc <- predict(cell$gc_fit, gc_wf, type = 'adjustedPeak') * cell$ploidy 
        cn <- NR_map_wcf / ER_gc * cell$ploidy # thus, in late S modal_CN windows will have cn == 4

        plot(1:length(cn), cn, pch = 19, cex = 0.4, col = col)
        abline(h = 0:6)
    })
}

fitCell <- function(cell_id, stage = "extract"){ 
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
