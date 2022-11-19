projectNameReactive <- function(sourceId){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        gsub('.scCNV.extract', '', getSourceFilePackageName(sourceId))
    })    
} 
normalizeDataReactive <- function(sourceId){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        dataFilePath <- getSourceFilePath(sourceId, "normalizeFile")
        startSpinner(session, message = "loading sample data")
        x <- readRDS(dataFilePath)
        setkey(x$colData, "cell_id")  
        x$qcPlotsDir <- expandSourceFilePath(sourceId, "qc_plots")
        if(!dir.exists(x$qcPlotsDir)) untar(
            getSourceFilePath(sourceId, "plotsArchive"), 
            exdir = x$qcPlotsDir
        )
        stopSpinner(session)
        x
    }) 
}

loadSampleRds <- function(sourceId){
    dataFilePath <- getSourceFilePath(sourceId, "normalizeFile")
    x <- readRDS(dataFilePath)
    setkey(x$colData, "cell_id")  
    x
}
loadSampleCommon <- function(cacheKey, keyObject, key, cacheObject, sourceId, ...){
    dprint("loadSampleCommon")
    x <- loadSampleRds(sourceId)
    x[c("rowRanges", "colData")]
}
loadSampleWorking <- function(cacheKey, keyObject, key, cacheObject, sourceId, chrom){
    dprint("loadSampleWorking")
    x <- loadSampleRds(sourceId)
    chrom_ <- chrom
    isAll <- chrom == "all"
    if(isAll){
        goodBins <- x$rowRanges[, bin_n]
    } else {
        goodBins <- x$rowRanges[chrom == chrom_, chrom_bin_n]
    }
    nGoodBins <- length(goodBins)
    allBins   <- 1:max(goodBins)
    badBins   <- allBins[!(allBins %in% goodBins)]
    nBadBins  <- length(badBins)    
    nCells <- nrow(x$colData)
    nullCell <- rep(NA,  nGoodBins) # we DO save all cells to allow dynamic cell rejection
    nullMatrix <- matrix(NA, nrow = nBadBins, ncol = nCells)
    order <- order(c(goodBins, badBins))
    sd <- sapply(x$cells, function(cell){
        i <- cell$hmm > 0
        sd((cell$cn[i] - cell$hmm[i]) / cell$hmm[i], na.rm = TRUE) # NA if NULL, etc.
    })
    getZ <- function(val) lapply(-1:1, function(increment){
        sapply(x$colData$cell_id, function(cell_id){
            cell <- x$cells[[cell_id]]
            if(is.null(cell[[val]])) return(nullCell)            
            w <- x$windows[[paste("w", cell$window_size, sep = "_")]]
            bins <- w[, if(isAll) TRUE else chrom == chrom_]
            x <- if(val == "hmm"){
                ifelse( # hmm as RGB intensities
                    cell[[val]][bins] == cell$modal_CN, 
                    0.5,
                    ifelse(
                        sign(cell[[val]][bins] - cell$modal_CN) == sign(increment),
                        1,
                        0
                    )
                )
            } else { # cn as CNC z-scores
                cn <- cell$modal_CN + increment
                (cell[[val]][bins] - cn) / cn / sd[cell_id]                
            }
            uncollapseVector(x, cell$window_size, nGoodBins)    
        })
    })
    padZ <- function(z) {
        if(nBadBins > 0) lapply(z, function(m) rbind(m, nullMatrix)[order, ]) else z
    }

    valueTypes <- c("cn", "hmm")
    windowMedians <- lapply(valueTypes, function(val){
        d <- sapply(x$colData$cell_id, function(cell_id){
            cell <- x$cells[[cell_id]]
            if(is.null(cell[[val]])) return(nullCell)            
            w <- x$windows[[paste("w", cell$window_size, sep = "_")]]
            bins <- w[, if(isAll) TRUE else chrom == chrom_]
            uncollapseVector(cell[[val]][bins], cell$window_size, nGoodBins)   
        })
        d <- apply(d, 1, median, na.rm = TRUE)
        lapply(x$windows, function(w) {
            n <- nGoodBins / w[chrom == chrom_, .N]
            collapseVector(d, n) / n
        })
    })
    names(windowMedians) <- valueTypes

    getValues <- function(val){
        x <- lapply(x$colData$cell_id, function(cell_id){
            cell <- x$cells[[cell_id]]
            if(is.null(cell[[val]])) return(NA)
            wx <- paste("w", cell$window_size, sep = "_")
            w <- x$windows[[wx]]
            bins <- w[, if(isAll) TRUE else chrom == chrom_]
            cell[[val]][bins] - windowMedians[[val]][[wx]]
        })
        names(x) <- x$colData$cell_id
        x
    }
    list(
        cell_id     = x$colData$cell_id,
        modal_CN    = x$colData$modal_CN,
        rejected    = x$colData$rejected,
        window_size = x$colData$window_size,
        windows = lapply(x$windows, function(w) w[chrom == chrom_, .(start, end)]),
        cn  = getValues("cn"),
        hmm = getValues("hmm"),
        z = list(
            cn  = padZ(getZ("cn")),
            hmm = padZ(getZ("hmm"))
        )
    )
}

#  $ 12 :List of 14
#   ..$ rejected   : logi FALSE
#   ..$ stage      : chr "extract"
#   ..$ pass       : num 1
#   ..$ modal_CN   : int 2
#   ..$ minBinCount: num 64
#   ..$ window_size: num 1
#   ..$ NR_map_w   : num [1:120949] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ gc_w       : num [1:120949] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ fit        :List of 8
#   .. ..$ nGcSteps     : num 100
#   .. ..$ gcIndexOffset: num -1
#   .. ..$ gcFractions  : num [1:60] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.0
# 9 ...
#   .. ..$ minGcIndex   : num 0
#   .. ..$ maxGcIndex   : num 59
#   .. ..$ theta        : num [1:60] 10 10 10 10 10 10 10 10 10 10 ...
#   .. ..$ mu           : num [1:60] 1.46 1.58 1.69 1.8 1.91 ...
#   .. ..$ peak         : int [1:60] 1 1 1 1 1 1 1 1 2 2 ...
#   .. ..- attr(*, "class")= chr "nbinomCountsGC2"
#   ..$ ER_gc      : num [1:120949] 2.92 2.92 2.92 2.92 2.92 ...
#   ..$ cn         : num [1:120949] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ hmm        : int [1:120949] 2 2 2 2 2 2 2 2 2 2 ...
#   ..$ percentile : num [1:120949] NA NA NA NA NA NA NA NA NA NA ...

# # expandWindows <- function(cell_id, data, key){ # go from window values to bin values, repeating window value over all bins
# #     cell <- cells[[cell_id]] # NOT data!
# #     w <- cell$window_size
# #     rollingRanges <- rollingRanges[[paste("w", w, sep = "_")]]
# #     chroms <- rollingRanges[reference_window == TRUE, chrom]
# #     unlist(rollingRanges[,
# #         sapply(data[[cell_id]][[key]][chroms == chrom[1]], function(x) rep(x, w))[1:.N],
# #         by = chrom
# #     ][, 2])
# # }
# # cn_w <- as.data.table(
# #     mclapply(constants$good_cell_ids, expandWindows, cells, "cn", mc.cores = env$N_CPU) 
# # )
