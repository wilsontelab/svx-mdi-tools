
    # CNV <- cell$windows$sequential$HMM != cell$ploidy

    # cell$windows$shape <- switch(
    #     env$SHAPE_CORRECTION,
    #     none = 1,
    #     batch = batchShapes[[as.character(cell$windowPower)]],
    #     cell = windows[, { # shape is relative to ER_modal_NA
    #         fI <- fitI[.I] # model is created from filtered windows
    #         nChromWindows <- sum(fI)
    #         if(nChromWindows < 10) rep(cell$windows$ER_modal_NA, .N) else { # insufficient points to fit
    #             isCnv <- CNV[.I][fI] # if <1/3 of chromosome is CNV, fit shape to the non-CNV portion only
    #             fracCNV <- sum(isCnv) / nChromWindows
    #             cnvI <- if(fracCNV < 1/3) !isCnv else TRUE
    #             i_f  <- i[fI][cnvI]
    #             i2_f <- i2[fI][cnvI]
    #             NR_wmhrf <- cell$windows$NR_wmhr[.I][fI][cnvI]
    #             fit <- lm(NR_wmhrf ~ i_f + i2_f)
    #             predict(fit, newdata = data.frame(i_f = i, i2_f = i2)) # final shape prediction occurs on ALL bins 
    #         }
    #     }, by = "chrom"][[2]] / cell$windows$ER_modal_NA
    # )
    # cell$windows$NR_wms <- with(cell$windows, { NR_wm / shape })
    # cell$windows$NR_wms <- with(cell$windows, { NR_wms * sum(NR_wm, na.rm = TRUE) / sum(NR_wms, na.rm = TRUE) }) # maintain total cell weights
    # plotWindows_counts(cell, "shapes", cell$windows$NR_wmhr[fitI],
    #                    windowI = fitI, shape = cell$windows$shape[fitI])
    
        # # re-solve the GC bias fit and squashed HMM using the reshaped NR values
    # cell$tmp$NR_wmshr <- with(cell, { windows$NR_wms * ploidy / (windows$sequential$HMM + windows$sequential$NAR) }) # thus, correct counts toward ploidy
    # col_w <- ifelse(fitI, defaultPointColor, rgb(1, 0, 0, 0.1))
    # cell <- fitGcBias(cell, "reshaped", cell$tmp$NR_wmshr, windowI = fitI, col_w = col_w)
    # plotWindows_counts(cell, "reshapedxx", cell$windows$NR_wms[fitI], windowI = fitI, col = "useRepColor")  

    # solveCompositeHmm(cell) # pass the baton
    
    # cell$repGcRatio <- NA
    # while(TRUE){
    #     cell <- setBestReplicationModel(cell, models, name)
    #     # solve an HMM to establish the (un)replicated genome spans for adjusting window NR prior to GC fit
    #     cell$windows$NAR <- if(cell$cellIsReplicating){
    #         emissProbs <- getRepEmissProbs(cell)
    #         hmm <- new_hmmEPTable(emissProbs, transProb = 1e-1, keys = windows$chrom[fitI])
    #         NAR <- keyedViterbi(hmm) - 1 # based on NR_wmhf
    #         round(NAR * cell$windows$HMM[fitI] / cell$ploidy, 0) # rescale NAR to hmm
    #     } else rep(0, sum(fitI))
    #     if(!cell$cellIsReplicating) break

    #     # check whether a putative replicating cell as the expected replication GC bias
    #     # if not, prevent the cell from being labeled as replicating
    #     windowIsReplicated <- cell$windows$NAR > 0        
    #     gc_replicated   <- windows[fitI][windowIsReplicated == TRUE,  mean(gc_fraction)]
    #     gc_unreplicated <- windows[fitI][windowIsReplicated == FALSE, mean(gc_fraction)]
    #     cell$repGcRatio <- gc_replicated / gc_unreplicated
    #     if(any(windowIsReplicated) && cell$repGcRatio > 1) break 
    #     models$allow <- models$modelType == "notReplicating"
    # }
    
    # use Q95 of window-to-window count deltas (NOT of raw window values) to determine 
# the window size that puts the majority of window counts between CN = ploidy +/- 1
# for most good cells, this default windowPower is > minWindowPower due to overdispersion
# normLagDiffQ_quantile  <- 0.5   # determined empirically
# normLagDiffQ_threshold <- 0.175 # determined empirically
# getNormLagDiffQ <- function(NR){ # 95%ile of the magnitude of the window-to-window count difference, normalized to local mean count
#     lagDiff <- abs(diff(NR)) 
#     lagMean <- (head(NR, -1) + tail(NR, -1)) / 2
#     quantile(lagDiff / lagMean, normLagDiffQ_quantile, na.rm = TRUE)
# }

# while(cell$windows$normLagDiffQ > normLagDiffQ_threshold && # stop when we get a nice, appropriately tight distribution

    # list( # one of these is chosen to become the working cell object
    #     NR_wm       = NR_wm,        
    #     ER_modal_NA = ER_modal_NA,
    #     ER_ploidy   = ER_modal_NA,
    #     RPA         = ER_modal_NA / env$PLOIDY, # later we may learn RPA is 2-fold off for replicating cells,
    #     normLagDiffQ = getNormLagDiffQ(NR_wmf)
    # )
    
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

    # cell$replicationModel$P_replicated_wf <- getP_replicated_wf(cell, gc_wf)
    # P_replicated_wf <- cell$replicationModel$P_replicated_wf # since cell is reset below
    # P_replicated_w <- getP_replicated_w(cell, gc_w) # a vector of reference GCI values

# getP_replicated_wf <- function(cell, gc_wf) with(cell, {
#     gci_wf <- round(gc_wf * nGcSteps, 0) # vector of filtered window GC indices
#     sapply(gcIndices, function(gci){ # all bins have at least somelikelihood of being called replicating
#         NAR <- windows$NAR[gci_wf == gci] # aggregate all windows with a given GC index
#         if(length(NAR) == 0) return(NA)
#         agg <- aggregate(NAR, list(NAR), length)
#         getN <- function(nar){
#             n <- agg[agg[[1]] == nar, 2]
#             if(length(n) == 0) 0 else n
#         }
#         pNAR <- sapply(0:ploidy, getN) / sum(agg[[2]]) # NAR's as observed for ploidy values in initial modeling
#         if(ploidy == 1){
#             sapply(0:4, function(cn) pmax(minPNar,
#                      if(cn == 0) 1    # one state, replication not meaningful
#                 else if(cn == 1) pNAR # CN==1, same as the NAR model
#                 else if(cn == 2) c(pNAR[1], minPNar, pNAR[2]) # CN=2,3,4, interpolate intermediate values
#                 else if(cn == 3) c(pNAR[1], minPNar, minPNar, pNAR[2])
#                 else             c(pNAR[1], minPNar, minPNar, minPNar, pNAR[2])
#             ))
#         } else { # ploidy == 2, typical
#             sapply(0:4, function(cn) pmax(minPNar,
#                      if(cn == 0) 1                       # one state, replication not meaningful
#                 else if(cn == 1) c(pNAR[1], 1 - pNAR[1]) # two states, unreplicated or replicated
#                 else if(cn == 2) pNAR                    # CN==2, same as the NAR model
#                 else if(cn == 3) c(pNAR[1], pNAR[2], pNAR[2], pNAR[3]) # CN=3,4
#                 else             c(pNAR[1], pNAR[2], pNAR[2], pNAR[2], pNAR[3])
#             ))
#         }
#     })
# })
# getP_replicated_w <- function(cell, gc_w)  with(cell$replicationModel, {
#     gci_w <- round(gc_w * nGcSteps, 0) # vector of all window GC indices
#     is <- 1:length(gci_w)    
#     nonNaGci <- which(!is.na(P_replicated_wf)) # GC indices with sufficient filtered windows to estimate P_unreplicated_wf
#     sapply(is, function(i){ # return the appropriate GCI per window (may have been adjusted from the actual GCI)
#         gci <- gci_w[i]
#         x <- P_replicated_wf[gci]
#         if(length(x) > 0 && !is.na(x)) return(gci)
#         delta <- abs(nonNaGci - gci)
#         nonNaGci[which(delta == min(delta))[1]]
#     })
# })
# commitReplicationProfile <- function(cell, gc_wf) with(cell, { # keep track of GC replication bias across cells
#     if(!cellIsReplicating) return(NULL)
#     gci_wf <- round(gc_wf * nGcSteps, 0) 
#     P_rep_fs_gc[[cell_id]] <<- list(
#         cell_id = cell_id,
#         windowPower = windowPower,
#         modelType = replicationModel$modelType,
#         fractionS = replicationModel$fractionS,
#         profile = sapply(gcIndices, function(gci){
#             NAR <- replicationModel$NAR[gci_wf == gci]
#             if(length(NAR) == 0) return(rep(NA, ploidy + 2))
#             agg <- aggregate(NAR, list(NAR), length)
#             N <- sum(agg[[2]])
#             c(N, sapply(0:ploidy, function(nar) {
#                 n <- agg[agg[[1]] == nar, 2]
#                 if(length(n) == 0) 0 else n
#             }) / N)
#         })
#     )  
# })
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
            # pNAR <- sapply(P_replicated_w, function(gci) P_replicated_wf[[gci]][[state$CN + 1]][state$NAR + 1])
            # dnbinom(NR_wms, mu = ER, size = theta) * pNAR
# P_rep_fs_gc <- list() # for collecting replicating cell profiles for later aggregation
