#----------------------------------------------------------------------
# server components for the normalizeBatch appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
normalizeBatchServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'normalizeBatch'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# sample summary table
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single",
    escape = FALSE,
    extraColumns = reactive({
        dt <- app$adjust$sampleSummaryColumns()
        dt[, .SD, .SDcols = c("cells", "keep", "reject", "override", "status")]
    })
)

#----------------------------------------------------------------------
# batch normalization actions
#----------------------------------------------------------------------
pendingSamples <- reactive({
    dt <- app$adjust$sampleSummaryColumns()
    req(dt)
    dt[batchFileExists == FALSE]  
})
output$executeUI <- renderUI({
    pendingSamples <- pendingSamples()
    bsButton(session$ns("execute"), "Run Pending Normalization(s)", 
             style = "primary", disabled = nrow(pendingSamples) == 0)
})
observeEvent(input$execute, {
    pendingSamples <- pendingSamples()
    for(i in seq_along(nrow(pendingSamples))){
        x <- loadSampleRds(pendingSamples[i, sourceId])
        # TODO: unlink old files ?
        for(chrom_ in x$constants$chroms[2 ]){
            goodBins <- x$rowRanges[chrom == chrom_, chrom_bin_n]
            nGoodBins <- length(goodBins)
            percentiles <- sapply(pendingSamples[i, cell_id][[1]], function(cell_id){
                cell <- x$cells[[cell_id]]
                wx <- paste("w", cell$window_size, sep = "_")
                w <- x$windows[[wx]]
                bins <- w[, chrom == chrom_]
                uncollapseVector(cell$percentile[bins], cell$window_size, nGoodBins)

# List of 14
#  $ rejected   : logi FALSE
#  $ stage      : chr "extract"
#  $ pass       : num 2
#  $ modal_CN   : int 2
#  $ minBinCount: num 82.3
#  $ window_size: num 17
#  $ NR_map_w   : num [1:7125] NA NA NA NA NA ...
#  $ gc_w       : num [1:7125] 0 0 0 0 0 ...
#  $ fit        :List of 8
#   ..$ nGcSteps     : num 100
#   ..$ gcIndexOffset: num -1
#   ..$ gcFractions  : num [1:56] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 .
# ..
#   ..$ minGcIndex   : num 0
#   ..$ maxGcIndex   : num 55
#   ..$ theta        : num [1:56] 10 10 10 10 10 10 10 10 10 10 ...
#   ..$ mu           : num [1:56] 47.9 48.1 48.3 48.5 48.7 ...
#   ..$ peak         : int [1:56] 43 43 43 43 43 44 44 44 44 44 ...
#   ..- attr(*, "class")= chr "nbinomCountsGC2"
#  $ ER_gc      : int [1:7125] 86 86 86 86 86 86 86 86 88 88 ...
#  $ cn         : num [1:7125] NA NA NA NA NA ...
#  $ hmm        : int [1:7125] 2 2 2 2 2 2 2 2 2 2 ...
#  $ percentile : num [1:7125] NA NA NA NA NA ...
#  $ theta      : num [1:7125] 10 10 10 10 10 ...

            })

            dstr(percentiles)




            # allBins   <- 1:max(goodBins)
            # badBins   <- allBins[!(allBins %in% goodBins)]
            # nBadBins  <- length(badBins)    
            # nCells <- nrow(x$colData)
            # nullCell <- rep(NA,  nGoodBins) # we DO save all cells to allow dynamic cell rejection
            # nullMatrix <- matrix(NA, nrow = nBadBins, ncol = nCells)
            # order <- order(c(goodBins, badBins))
            # sd <- sapply(x$cells, function(cell){
            #     i <- cell$hmm > 0
            #     sd((cell$cn[i] - cell$hmm[i]) / cell$hmm[i], na.rm = TRUE) # NA if NULL, etc.
            # })
        }

#  $ bin_size          : int 20000
#  $ chroms            : chr [1:21(1d)] "chr1" "chr2" "chr3" "chr4" ...
#  $ genome_bins       : int 136288
#  $ num_bins_per_chrom: int [1:21(1d)] 9774 9106 8002 7826 7592 7487 7273 6471 62
# 30 6535 ...
#  $ num_cells         : int 394
#  $ num_chroms        : int 21
#  $ num_nodes         : int 787

# [1] "env"       "constants" "metadata"  "rowRanges" "windows"   "colData"       
# [7] "cells"



        # getZ <- function(val) lapply(-1:1, function(increment){
        #     sapply(x$colData$cell_id, function(cell_id){
        #         cell <- x$cells[[cell_id]]
        #         if(is.null(cell[[val]])) return(nullCell)            
        #         w <- x$windows[[paste("w", cell$window_size, sep = "_")]]
        #         bins <- w[, if(isAll) TRUE else chrom == chrom_]
        #         x <- if(val == "hmm"){
        #             ifelse( # hmm as RGB intensities
        #                 cell[[val]][bins] == cell$modal_CN, 
        #                 0.5,
        #                 ifelse(
        #                     sign(cell[[val]][bins] - cell$modal_CN) == sign(increment),
        #                     1,
        #                     0
        #                 )
        #             )
        #         } else { # cn as CNC z-scores
        #             cn <- cell$modal_CN + increment
        #             (cell[[val]][bins] - cn) / cn / sd[cell_id]                
        #         }
        #         uncollapseVector(x, cell$window_size, nGoodBins)    
        #     })
        # })
        # padZ <- function(z) {
        #     if(nBadBins > 0) lapply(z, function(m) rbind(m, nullMatrix)[order, ]) else z
        # }

        # valueTypes <- c("cn", "hmm")
        # windowMedians <- lapply(valueTypes, function(val){
        #     d <- sapply(x$colData$cell_id, function(cell_id){
        #         cell <- x$cells[[cell_id]]
        #         if(is.null(cell[[val]])) return(nullCell)            
        #         w <- x$windows[[paste("w", cell$window_size, sep = "_")]]
        #         bins <- w[, if(isAll) TRUE else chrom == chrom_]
        #         uncollapseVector(cell[[val]][bins], cell$window_size, nGoodBins)   
        #     })
        #     d <- apply(d, 1, median, na.rm = TRUE)
        #     lapply(x$windows, function(w) {
        #         n <- nGoodBins / w[chrom == chrom_, .N]
        #         collapseVector(d, n) / n
        #     })
        # })
        # names(windowMedians) <- valueTypes

        # getValues <- function(val){
        #     x <- lapply(x$colData$cell_id, function(cell_id){
        #         cell <- x$cells[[cell_id]]
        #         if(is.null(cell[[val]])) return(NA)
        #         wx <- paste("w", cell$window_size, sep = "_")
        #         w <- x$windows[[wx]]
        #         bins <- w[, if(isAll) TRUE else chrom == chrom_]
        #         cell[[val]][bins] # - windowMedians[[val]][[wx]]
        #     })
        #     names(x) <- x$colData$cell_id
        #     x
        # }
        # list(
        #     cell_id     = x$colData$cell_id,
        #     modal_CN    = x$colData$modal_CN,
        #     rejected    = x$colData$rejected,
        #     window_size = x$colData$window_size,
        #     windows = lapply(x$windows, function(w) w[chrom == chrom_, .(start, end)]),
        #     windowMedians = windowMedians,
        #     cn  = getValues("cn"),
        #     hmm = getValues("hmm"),
        #     z = list(
        #         cn  = padZ(getZ("cn")),
        #         hmm = padZ(getZ("hmm"))
        #     )
        # )


        # dstr(x)
    }

})

# userOverrides[[sourceId]][[input$cellKeepReject$cell_id]]
# dstr(app$adjust$outcomes$userOverrides())


#----------------------------------------------------------------------
# aggregate plots of all good cells
#----------------------------------------------------------------------
sampleData <- reactive({
    sid <- sourceId()
    req(sid)
    d <- loadSampleRds(sid)
    d$good_cell_ids <- pendingSamples()[sourceId == sid, cell_id][[1]]
    d$colData <- d$colData[d$good_cell_ids]
    d$constants$min_window_size <- min(d$colData$window_size, na.rm = TRUE)
    d
})
windowSize <- staticPlotBoxServer(
    "windowSize",
    #----------------------------
    maxHeight = "400px",
    immediate = TRUE,
    title = FALSE,
    #----------------------------
    envir = parent.frame(),
    settings = NULL,
    create = function(){

        # req(FALSE)
        
        d <- sampleData()
        d <- aggregate(d$colData$window_size, list(d$colData$window_size), length)        
        par(mar = c(4.1, 4.1, 0.1, 0.1))
        windowSize$initializeFrame(
            # title = pendingSamples[sourceId == sid, Project],
            xlim = range(d[[1]], na.rm = TRUE),
            ylim = range(d[[2]], na.rm = TRUE),
            xlab = "Window Size (# of 20kb bins)",
            ylab = "# of Cells"
        )
        windowSize$addPoints(
            x = d[[1]],
            y = d[[2]],
            typ = "h"
        )
    }
)

expandWindows <- function(d, val){
    dmsg(paste("expanding", val))
    do.call(rbind, lapply(d$constants$chroms, function(chrom_){
        nGoodBins <- d$rowRanges[chrom == chrom_, .N]
        x <- sapply(d$good_cell_ids, function(cell_id){
            cell <- d$cells[[cell_id]]
            wx <- paste("w", cell$window_size, sep = "_")
            bins <- d$windows[[wx]][, chrom == chrom_]  
            workingWindowSize <- cell$window_size / d$constants$min_window_size
            workingNBins <- ceiling(nGoodBins / d$constants$min_window_size)
            if(workingWindowSize == 1) cell[[val]][bins]
            else uncollapseVector(cell[[val]][bins], workingWindowSize, workingNBins)
        })
        dprint(paste(chrom_, nrow(x)))
        x
    }))    
}
getRows_nbinomCountsGC2 <- function(nb, fractionGC){
    pmax(nb$minGcIndex, pmin(nb$maxGcIndex, round(fractionGC * nb$nGcSteps, 0) - nb$gcIndexOffset))
}
predict_nbinomCountsGC2 <- function(nb, fractionGC){
    nb$mu[getRows_nbinomCountsGC2(nb, fractionGC)]
}

batchEffect <- staticPlotBoxServer(
    "batchEffect",
    #----------------------------
    maxHeight = "400px",
    immediate = FALSE,
    title = FALSE,
    #----------------------------
    envir = parent.frame(),
    settings = NULL,
    create = function(){
        d <- sampleData()

        sid <- sourceId()
        sidG <- tryCatch({get("sourceIdG", envir = .GlobalEnv)}, error = function(e) "XXX")
        forceAll <- sid != sidG
        assign("sourceIdG", sid, envir = .GlobalEnv)
        
        if(forceAll || !exists("hmmG", envir = .GlobalEnv)) assign("hmmG", round(expandWindows(d, "hmm"), 0), envir = .GlobalEnv)
        if(forceAll || !exists("medianHmmG", envir = .GlobalEnv)) assign("medianHmmG", round(apply(hmmG, 1, median, na.rm = TRUE), 0), envir = .GlobalEnv)
        workingCells <- hmmG == medianHmmG

        workingWindowSizes <- sort(unique(d$colData$window_size))
        windowSizePairs <- as.data.table(expand.grid(workingWindowSizes, workingWindowSizes))
        setnames(windowSizePairs, c("ws1", "ws2"))
        windowSizePairs$key <- apply(windowSizePairs, 1, paste, collapse = ":")         
        windowSizePairs[, weight := pmin(ws1, ws2) / pmax(ws1, ws2)]        
        setkey(windowSizePairs, "key")

        # if(!exists("percentileG", envir = .GlobalEnv)) assign("percentileG", expandWindows(d, "percentile"), envir = .GlobalEnv)
        if(forceAll || !exists("cnG", envir = .GlobalEnv)) assign("cnG", expandWindows(d, "cn"), envir = .GlobalEnv)
        ws2 <- d$colData[d$good_cell_ids, window_size]

        workingWindowSizes <- 2

        # if(!exists("meanPercentileG", envir = .GlobalEnv)) assign("meanPercentileG", lapply(workingWindowSizes, function(ws1){
        #     x <- sapply(1:nrow(percentileG), function(i){

        #         j <- workingCells[i, ] & ws2 == ws1

        #         median(percentileG[i, j], na.rm = TRUE)

        #         # keys <- paste(ws1, ws2[j], sep = ":")
        #         # weighted.mean(percentileG[i, j], windowSizePairs[keys, weight], na.rm = TRUE)
        #     })
        #     wws1 <- ws1 / d$constants$min_window_size
        #     if(wws1 == 1) x else collapseVector(x, wws1) / wws1
        # }), envir = .GlobalEnv)
        # names(meanPercentileG) <- as.character(workingWindowSizes) 

        if(forceAll || !exists("meanCncG", envir = .GlobalEnv)) assign("meanCncG", {
            dmsg("meanCncG")
            lapply(workingWindowSizes, function(ws1){
                x <- sapply(1:nrow(cnG), function(i){
                    j <- workingCells[i, ]
                    if(sum(j) == 0) return(NA)
                    cn <- cnG[i, j]
                    if(all(is.na(cn))) return(NA)
                    keys <- paste(ws1, ws2[j], sep = ":") # cannot embed directly in windowSizePairs[]
                    weighted.mean(cn - medianHmmG[i], windowSizePairs[keys, weight], na.rm = TRUE)
                })
                wws1 <- ws1 / d$constants$min_window_size
                if(wws1 == 1) x else collapseVector(x, wws1) / wws1
            })
        }, envir = .GlobalEnv)
        names(meanCncG) <- as.character(workingWindowSizes) 

        dmsg("plotting")

        # cell_id <- "7" #"2" "234" "91" "288"
        cell_id <-d$colData[d$good_cell_ids][window_size == workingWindowSizes, sample(cell_id, 1)]
        cellI <- which(d$good_cell_ids == cell_id)
        cell <- d$cells[[cell_id]]
        ws1 <- cell$window_size 
        dprint(cell_id)

        # dd <- round(meanPercentileG[[as.character(ws1)]], 2)
        # dd <- aggregate(dd, list(dd), length)

        # par(mar = c(4.1, 4.1, 0.1, 0.1))
        # windowSize$initializeFrame(
        #     # title = pendingSamples[sourceId == sid, Project],
        #     xlim = range(dd[[1]]),
        #     ylim = c(0, max(dd[[2]])),
        #     xlab = "meanPercentileG",
        #     ylab = "Frequency"
        # )
        # windowSize$addLines(
        #     x = dd[[1]],
        #     y = dd[[2]],
        #     col = "blue"
        # )
        # abline(v = 0.5)

        # meanPrc <- meanPercentileG[[as.character(ws1)]]
        # cellPrc <- percentileG[, cellI]
        # wws1 <- ws1 / d$constants$min_window_size
        # if(wws1 != 1) cellPrc <- collapseVector(cellPrc, wws1) / wws1
        # fit <- lm(cellPrc ~ meanPrc)
        # fitPrc <- predict(fit)
        # rpa <- predict_nbinomCountsGC2(cell$fit, cell$gc_w)
        # refPrc <- median(fitPrc, na.rm = TRUE)
        # excess <- qnbinom(fitPrc, size = cell$theta, mu = rpa) -
        #           qnbinom(refPrc, size = cell$theta, mu = rpa)
        # ER1 <- rpa * cell$modal_CN
        # ER2 <- (rpa + excess) * cell$modal_CN
        # cn1 <- cell$NR_map_w / ER1 * cell$modal_CN
        # cn2 <- cell$NR_map_w / ER2 * cell$modal_CN

        meanCnc <- meanCncG[[as.character(ws1)]]
        cnc <- cnG[, cellI] - medianHmmG
        modalBins <- workingCells[, cellI]        
        wws1 <- ws1 / d$constants$min_window_size
        if(wws1 != 1) {
            cnc <- collapseVector(cnc, wws1) / wws1
            modalBins <- collapseVector(modalBins, wws1) > 0
        }
        fitBins <- which(meanCnc > -1 & meanCnc < 1 & modalBins)
        fit <- lm(cnc ~ meanCnc, data = data.frame(cnc = cnc[fitBins], meanCnc = meanCnc[fitBins]))
        cncFit <- cnc - predict(fit, newdata = data.frame(meanCnc == meanCnc))

        par(mar = c(4.1, 4.1, 0.1, 0.1))
        bin <- 1:length(cnc)
        windowSize$initializeFrame(
            # title = pendingSamples[sourceId == sid, Project],
            xlim = range(bin),
            ylim = c(-2, 2),
            xlab = "bin",
            ylab = "CNC"
        )
        x <- c(bin, bin)
        # y <- c(cnc, cncFit)
        y <- c(rep(NA, length(cnc)), cncFit)
        col <- c(rep("blue", length(cnc)), rep("red", length(cncFit)))
        N <- length(x)
        order <- sample(sample(N), min(N, 100000))
        windowSize$addPoints(
            x = x[order],
            y = y[order],
            pch = 16,
            cex = 0.25,
            col = col[order]
        )

        # par(mar = c(4.1, 4.1, 0.1, 0.1))
        # cnc    <- round(cnc,    1)
        # cncFit <- round(cncFit, 1)
        # cnc    <- aggregate(cnc,    list(cnc),    length)
        # cncFit <- aggregate(cncFit, list(cncFit), length)        
        # windowSize$initializeFrame(
        #     # title = pendingSamples[sourceId == sid, Project],
        #     xlim = c(-2, 2),
        #     ylim = c(0, max(cnc[[2]], cncFit[[2]])),
        #     xlab = "CNC",
        #     ylab = "Frequency"
        # )
        # windowSize$addLines(
        #     x = cnc[[1]],
        #     y = cnc[[2]],
        #     col = "blue"
        # )
        # windowSize$addLines(
        #     x = cncFit[[1]],
        #     y = cncFit[[2]],
        #     col = "red"
        # )

        # par(mar = c(4.1, 4.1, 0.1, 0.1))
        # windowSize$initializeFrame(
        #     # title = pendingSamples[sourceId == sid, Project],
        #     xlim = c(-1, 1),
        #     ylim = c(-2, 2),
        #     xlab = "Mean CNC",
        #     ylab = "CNC"
        # )
        # windowSize$addPoints(
        #     x = meanCnc,
        #     y = cnc,
        #     pch = 16, 
        #     cex = 0.25
        # )
        # coef <- coef(fit)
        # abline(coef[1], coef[2])
        # windowSize$addPoints(
        #     x = meanCnc,
        #     y = cncFit,
        #     pch = 16, 
        #     cex = 0.25,
        #     col = "red"
        # )
        # fit <- lm(cnc ~ meanCnc, data = data.frame(cnc = cncFit[fitBins], meanCnc = meanCnc[fitBins]))
        # coef <- coef(fit)
        # abline(coef[1], coef[2], col = "red")


        # par(mar = c(4.1, 4.1, 0.1, 0.1))
        # x <- meanPrc
        # y <- percentileG[, cellI]
        # wws1 <- ws1 / d$constants$min_window_size
        # if(wws1 != 1) y <- collapseVector(y, wws1) / wws1
        # windowSize$initializeFrame(
        #     # title = pendingSamples[sourceId == sid, Project],
        #     xlim = c(0.35, 0.7),
        #     ylim = range(y, na.rm = TRUE),
        #     xlab = "Mean Percentile",
        #     ylab = "Percentile"
        # )
        # windowSize$addPoints(
        #     x = x,
        #     y = y,
        #     pch = 16, 
        #     cex = 0.25
        # )
        # i <- which(x > 0.4 & x < 0.65)
        # fit <- lm(y ~ x, data = data.frame(x =x[i], y = y[i]))
        # coef <- coef(fit)
        # abline(coef[1], coef[2])
        # abline(0, 1, col = "blue")
        # y <- y - (predict(fit, newdata = data.frame(x = x, y = y)) - 0.5)
        # # windowSize$addPoints(
        # #     x = x,
        # #     y = y,
        # #     pch = 16, 
        # #     cex = 0.25,
        # #     col = "red"
        # # )
        # fit <- lm(y ~ x, data = data.frame(x =x[i], y = y[i]))
        # coef <- coef(fit)
        # abline(coef[1], coef[2], col = "red")
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
