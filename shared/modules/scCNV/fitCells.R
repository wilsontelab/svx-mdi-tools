# determine the least number of reads per window for a cell's data to support a robust HMM
# places ~96% of windows within (modal_CN +/- 0.5) * replicationFactor==[1,2]
# sapply(1:4, getMinWindowCount) => 16 64 144 256; thus, usually 64 reads/window for interphase diploid cells
getMinWindowCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
    (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# set the per-cell minimum window size as the number of bins needed to obtain a mean raw count >=minWindowCount
# cells with too-low counts will still be analyzed at MAX_WINDOW_POWER but will likely fail QC in getCellDefaults
windowSizes <- 2 ** windowPowers
bins <- rowRanges$autosome & rowRanges$mappability >= env$MIN_MAPPABILITY
getMinWindowPower <- function(cell_id, minWindowCount){
    NR_raw_b <- raw_counts[[cell_id]][bins]
    NR_avg_b <- mean(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_avg_b) || NR_avg_b == 0) return(env$MAX_WINDOW_POWER)
    windowSize <- ceiling(minWindowCount / NR_avg_b) # in number of bins, not bp
    if(is.na(windowSize) || windowSize > 2 ** env$MAX_WINDOW_POWER) return(env$MAX_WINDOW_POWER)
    windowSize <- windowSizes[min(which(windowSizes >= windowSize))] # thus, allow 2**(0:MAX_WINDOW_POWER) bins per window
    log2(windowSize)
}

# use sd of window-to-window count deltas (NOT sd of raw window values) to determine 
# the window size that puts the majority of window counts between CN = ploidy +/- 1
# flag cells that cannot achieve this metric as likely unusable cells
# for most good cells, windowPower is > minWindowPower due to overdispersion
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
getSdLagDiff <- function(cell_id, windowPower){
    windowSize <- 2 ** windowPower
    windows <- windows[[windowPower + 1]]
    mappability <- windows[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_map_w <- unname(unlist(sapply(constants$chrom, function(chrom){
        collapseVector(
            raw_counts[[cell_id]][rowRanges$chrom == chrom], 
            windowSize
        ) / mappability[windows$chrom == chrom]
    })))
    NR_map_wf <- excludeOutliers(NR_map_w, min = 1) # thus, chrom-to-chrom, baseline-to-CNV, and rare bad windows won't have undue influence
    ER_modal_CN <- peakValue(NR_map_wf)    # ER_modal_CN may NOT be the same as ER_ploidy
    readsPerAllele <- ER_modal_CN / env$PLOIDY # later we may learn readsPerAllele is 2-fold off for replicating cells 
    lagDiff <- excludeOutliers(diff(NR_map_wf), low = 0, high = 0.95) # further minimize the impact of segment transitions
    list( # one of these is chosen to become the initial cellDefaults object
        cell_id = cell_id,
        windowPower = windowPower,
        replicating = FALSE,  
        modal_CN = env$PLOIDY, # subject to change if late S-phase
        ploidy = env$PLOIDY,
        ER_modal_CN = ER_modal_CN,
        ER_ploidy = ER_modal_CN,
        readsPerAllele = readsPerAllele,
        sdLagDiff = sd(lagDiff) / readsPerAllele * env$PLOIDY,
        NR_map_w = NR_map_w
    )
}
getCellDefaults <- function(cell_id, minWindowPower){
    windowPower <- minWindowPower
    x <- getSdLagDiff(cell_id, windowPower)
    while(x$sdLagDiff > 1 && windowPower < env$MAX_WINDOW_POWER){
        windowPower <- windowPower + 1
        x <- getSdLagDiff(cell_id, windowPower)
    }
    x$keep <- x$sdLagDiff <= 1
    x$minWindowPower <- minWindowPower
    x
}

# working only at the optimized windowPower for the given cell, normalize it's data for aneuploidy and replication
getReplicationModel <- function(cell, peakIsReplicated){
    ER_replicated   <- cell$ER_modal_CN * (if(peakIsReplicated) 1 else 2)
    ER_unreplicated <- ER_replicated / 2
    sd_replicated   <- (ER_replicated - ER_unreplicated) / (if(peakIsReplicated) 2 else 4) # since window fitting above assumed modal_CN == ploidy
    sd_unreplicated <- sd_replicated / sqrt(2) # since sampling variance is proportional to mean
    maxNR <- ER_replicated   + 2 * sd_replicated
    minNR <- ER_unreplicated - 2 * sd_unreplicated
    list(
        ER_unreplicated = ER_unreplicated,
        sd_unreplicated = sd_unreplicated,
        ER_replicated   = ER_replicated,
        sd_replicated   = sd_replicated,
        minNR = minNR,
        maxNR = maxNR,
        model = data.table(
            P_unreplicated = dnorm(cell$NR_map_wc, mean = ER_unreplicated, sd = sd_unreplicated), # TODO: adjust the probability of being early replicating based on GC
            P_replicated   = dnorm(cell$NR_map_wc, mean = ER_replicated,   sd = sd_replicated)  # here, i.e., dnorm() * pEarly_gc, etc. (nearly just fracGC?)
        )
    )    
}
adjustReplicationModel <- function(model, fractionS){ # update model rep/unrep bin emisssion probabilities for fractionS in preparation for HMM
    model[, ":="(
        P_unreplicated = P_unreplicated * (1 - fractionS),        
        P_replicated   = P_replicated   * fractionS
    )]
    model
}
getModelLogLikelihood <- function(P_replicated, P_unreplicated, fractionS){
    x <- log(P_unreplicated * (1 - fractionS) + P_replicated * fractionS)
    x[x == -Inf] <- -800 # suppress outliers      
    sum(x)
}
applyCnReplicationCorrections <- function(cell){ # cellDefaults called just cell here for readability
    windows <- windows[[cell$windowPower + 1]]
    chroms  <- windows[, unique(chrom)]

    # determine a first crude aneuploidy copy number estimate for each whole chromosome
    maxModelCN <- 5
    cns <- 0:maxModelCN # CN==5 is not trustworthy, the true value could be (much) higher
    NR_map_wf <- excludeOutliers(cell$NR_map_w, min = 2)
    RPA_wf <- NR_map_wf / cell$ploidy
    NR_wf_cn <- lapply(cns, "*", RPA_wf)
    cell$chromCn <- lapply(chroms, function(chrom){
        NR_chrom <- excludeOutliers(cell$NR_map_w[windows$chrom == chrom], min = 2)
        if(length(NR_chrom) <= 1 || median(NR_chrom) <= 1) return(0)
        p <- sapply(NR_wf_cn, function(NR_g) wilcox.test(NR_g, NR_chrom, exact = FALSE)$p.value)
        cns[which.max(p)] # the copy number that gives the best distribution match between chrom and genome
    })
    names(cell$chromCn) <- chroms
    cell$chromCnIsTrustworthy <- cell$chromCn[windows$chrom] %between% c(1, 4)

    # use aneuploidy-corrected counts to establish replication models consistent with peak +/- replication offset from it
    cell$NR_map_wc <- unname(unlist(sapply(chroms, function(chrom){ # rescale aneuploid chromosomes (including chrX) to ploidy
        cell$NR_map_w[windows$chrom == chrom] * cell$ploidy / cell$chromCn[[chrom]]
    })))
    peakIsUnreplicated_model <- getReplicationModel(cell, FALSE) # first model handles unreplicated to early stage replication where peak CN typically == ploidy
    peakIsReplicated_model   <- getReplicationModel(cell, TRUE)  # second model handles late stage replication where peak CN typically == 2 * ploidy
    repFitI <- with(cell, { 
        chromCnIsTrustworthy & 
        !is.na(NR_map_wc) & 
        NR_map_wc >= peakIsReplicated_model$minNR & 
        NR_map_wc <= peakIsUnreplicated_model$maxNR 
    })
    NR_map_wcf <- cell$NR_map_wc[repFitI]
    replicationModels <- rbind(
        do.call(rbind, lapply(c(0, seq(0.05, 0.75, 0.01)), function(fractionS){ # explore replication weighting over 1% increments of S-phase
            data.table(                                                         # allow each possible peak location in mid-S
                peakIsReplicated = FALSE,
                fractionS = fractionS,
                cellIsReplicating = fractionS > 0,
                logLikelihood = peakIsUnreplicated_model$model[repFitI, getModelLogLikelihood(P_replicated, P_unreplicated, fractionS)]          
            )
        })),
        do.call(rbind, lapply(seq(0.25, 0.95, 0.01), function(fractionS){ # don't use 100% replicated, call that 0%
            data.table(
                peakIsReplicated = TRUE,
                fractionS = fractionS,
                cellIsReplicating = TRUE,
                logLikelihood = peakIsReplicated_model$model[repFitI, getModelLogLikelihood(P_replicated, P_unreplicated, fractionS)]          
            )
        }))
    )

    # select the best replication model, a combination of peak location (rep or unrep) and fractionS
    cell$replicationModel <- as.list(replicationModels[logLikelihood == max(logLikelihood)][fractionS == min(fractionS)][1])
    cell$modal_CN       <- cell$ploidy * (if(cell$replicationModel$peakIsReplicated) 2 else 1)
    cell$ER_ploidy <- cell$ER_modal_CN / (if(cell$replicationModel$peakIsReplicated) 2 else 1)
    cell$readsPerAllele <- cell$ER_modal_CN / cell$modal_CN # all not yet subjected to gc fit   
    cell$replicationModel$model <- with(cell$replicationModel, { adjustReplicationModel(
        if(peakIsReplicated) peakIsReplicated_model$model else peakIsUnreplicated_model$model,
        fractionS
    )})

    # solve an HMM to establish the (un)replicated genome spans for adjusting window NR prior to GC git
    emissProbs <- log(as.matrix(cell$replicationModel$model[repFitI]))
    hmm <- new_hmmEPTable(emissProbs, transProb = 1e-1, keys = windows$chrom[repFitI])
    cell$replicationModel$windowIsReplicated <- rep(NA, length(repFitI))
    cell$replicationModel$windowIsReplicated[repFitI] <- as.logical(keyedViterbi(hmm) - 1)
    # TMP_DIR <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV"
    # pngFile <- file.path(TMP_DIR, paste0("GC_bias", ".png"))
    # png(pngFile, width = 6, height = 3, units = "in", pointsize = 8, res = 300, type = "cairo")
    # dev.off()
    # stop("XXXXXXXXXXXXXXX")
    cell
}

# fit a GC bias correction to the normalized window data and solve for CNVs by HMM
fitCellGcHmm <- function(cell){
    windows <- windows[[cell$windowPower + 1]]
    repFitI <- !is.na(cell$replicationModel$windowIsReplicated)
    NR_map_wcr <- with(cell, { NR_map_wc / ifelse(replicationModel$windowIsReplicated, 2, 1) }) # thus, these are adjusted toward ploidy
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
    minWindowCount <- getMinWindowCount(env$PLOIDY)
    minWindowPower <- getMinWindowPower(cell_id, minWindowCount)
    cell <- getCellDefaults(cell_id, minWindowPower)

    # TODO: run this twice
    # first time, use crude CN correction
    # second time, pass the CN obtained from HMM in first round
    # also pass the GC fit from 1st round to the 2nd round?
    # consider also resetting window power using GC fit 2nd time?

    # espeically if windowPower setting is unreliable, consider doing this for several windowPowers per cell? min:default

    cell <- applyCnReplicationCorrections(cell)
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
    plot(cell$fit, cell$gc_w, cell$NR_map_w, cell$modal_CN, !cell$keep)

    # plot NR_map_w vs. window index, i.e., pre-normalization
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plotCellByWindow(cell$NR_map_w, "# Reads", ylim = c(0, cell$ER_ploidy * 3))

    # plot CN vs. window index, i.e., post-normalization
    plotCellByWindow(cell$cn, "CN", ylim = c(0, 6))
    abline(h = 0:4, col = "grey")
    if(!is.null(cell$hmm)) lines(1:length(cell$hmm), cell$hmm, col = "red")

    dev.off()
}

# fitCell_ <- function(cell_id, modal_CN, minBinCount, windowPower, pass, stage = "extract"){ 

#     # calculate window sums and correct for mappability
#     window_size <- 2 ** windowPower
#     windows <- windows[[paste("w", window_size, sep = "_")]]   
#     mappability <- windows[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
#     NR_map_w <- unname(unlist(sapply(constants$chrom, function(chrom){
#         collapseVector(
#             raw_counts[[cell_id]][rowRanges$chrom == chrom], 
#             window_size
#         ) / mappability[windows$chrom == chrom]
#     })))

#     # process required bin and cell information
#     gc_w  <- windows$gc_fraction
#     gc_wa <- gc_w[windows$autosome]
#     NR_map_wa <- NR_map_w[windows$autosome]
#     autosomes <- windows$chrom[windows$autosome]

#     # perform an initial GC bias fit that assumes the same CN for all autosomes
#     fit <- new_nbinomCountsGC2(NR_map_wa, gc_wa, binCN = modal_CN)

#     # use the initial fit to solve an initial CN estimate for all autosomes
#     hmm <- viterbi(fit, NR_map_wa, gc_wa, asRle = FALSE,
#                    chroms = autosomes, transProb = env$TRANSITION_PROBABILITY)
#     cn_estimate <- ifelse(hmm$cn == hmm$maxCN | hmm$cn == 0, NA, hmm$cn) # bins at maxCN are unreliable as they might be >maxCN

#     # revise to a final GC bias fit using the initial copy number estimates
#     fit <- new_nbinomCountsGC2(NR_map_wa, gc_wa, binCN = cn_estimate)

#     # solve a final CN estimate for all chromosomes (not just autosomes)
#     ER_gc <- predict(fit, gc_w, type = 'adjustedPeak') * modal_CN # use peak for visualization, unless it is unreliable   
#     cn <- NR_map_w / ER_gc * modal_CN 
#     hmm <- viterbi(fit, NR_map_w, gc_w, asRle = FALSE, 
#                    chroms = windows$chrom, transProb = env$TRANSITION_PROBABILITY)

#     # return our results
#     list(
#         rejected = TRUE, # caller must set to FALSE if keeping the cell
#         stage = stage,
#         pass = pass,
#         modal_CN = modal_CN,
#         minBinCount = minBinCount,
#         window_size = window_size,
#         NR_map_w = NR_map_w,            
#         gc_w = gc_w,
#         fit = fit, 
#         ER_gc = ER_gc,  
#         cn = cn,              
#         hmm = hmm$cn,
#         percentile = cumprob(fit, NR_map_w, gc_w, binCN = hmm$cn),         
#         theta = predict(fit, gc_w, type = 'theta')
#     )
# }

    # # step up to larger and larger windows as needed to account for overdispersion
    # # expanding window sizes by a factor of two leads to predictable cell-to-cell bin relationships
    # pass <- 0
    # windowPower <- minWindowPower - 1 
    # while(x$rejected && windowPower < env$MAX_WINDOW_POWER){
    #     pass <- pass + 1
    #     windowPower <- windowPower + 1
    #     x <- fitCell_(cell_id, modal_CN, minBinCount, windowPower, pass, stage)
    #     dcn <- x$cn[x$hmm == x$modal_CN] - x$modal_CN
    #     sdDcn <- sd(dcn, na.rm = TRUE)
    #     x$rejected <- is.na(sdDcn) || sdDcn * env$N_SD_HALFCN > 1 # deliberately 2-fold more permissive than the stated target
    # }

    # # return our result
    # x




    # # establish a density power for creating the initial CN HMM models, to ensure a reasonably smooth density
    # NR_map_wcf_int <- as.integer(round(NR_map_wcf, 0))
    # NR_model <- NR_map_wcf_int
    # nDensityValues <- diff(range(NR_model)) + 1
    # densityPower <- 0
    # while(nDensityValues > 100){
    #     densityPower <- densityPower + 1
    #     NR_model <- round(NR_map_wcf_int / 2**densityPower, 0)
    #     nDensityValues <- diff(range(NR_model)) + 1
    # }
    # maxPossibleModelValue <- round(max(NR_map_wcf) * maxModelCN / cell$ploidy / 2**densityPower, 0)
    # nullP <- rep(1e-6, maxPossibleModelValue + 1) # add 1 to account for NR=zero
    # NR_density <- round(cell$NR_map_w / 2**densityPower, 0)

    # # use the density power to create emmission probabilities for an initial HMM 
    # emissProbs <- sapply(cns, function(cn){
    #     if(cn == 0) cn <- 0.05 # allow for occasional stray reads at CN==0
    #     NR_model <- round(NR_map_wcf * cn / cell$ploidy / 2**densityPower, 0)
    #     agg <- aggregate(NR_model, list(NR_model), length)
    #     p <- nullP
    #     p[agg[[1]] + 1] <- agg[[2]] / sum(agg[[2]])
    #     log(p[NR_density + 1])
    # })

    # # solve the HMM and use it to refine the model to be used for replication assessment
    # hmm <- new_hmmEPTable(emissProbs, transProb = 1e-6, keys = windows$chrom)
    # cell$hmm1 <- keyedViterbi(hmm) - 1

    # TMP_DIR <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/shared/modules/scCNV"
    # pngFile <- file.path(TMP_DIR, paste0("hmm1", ".png"))
    # png(pngFile, width = 6, height = 3, units = "in", pointsize = 8, res = 300, type = "cairo")
    # plot(1:length(cell$NR_map_w), cell$NR_map_w, pch = 16, cex = 0.25)
    # lines(1:length(cell$NR_map_w), cell$hmm1 * cell$readsPerAllele, col = "red")
    # dev.off()

