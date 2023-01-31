# normalizeDataReactive <- function(sourceId, ...){
#     reactive({
#         sourceId <- sourceId()
#         req(sourceId)
#         normalizeFilePath <- getSourceFilePath(sourceId, "normalizeFile")
#         startSpinner(session, message = "loading sample data")
#         x <- readRDS(normalizeFilePath)
#         setkey(x$colData, "cell_id")  
#         x$qcPlotsDir <- expandSourceFilePath(sourceId, "qc_plots")
#         if(!dir.exists(x$qcPlotsDir)) dir.create(x$qcPlotsDir)
#         stopSpinner(session)
#         x
#     }) 
# }
# loadSampleRds <- function(sourceId){
#     normalizeFilePath <- getSourceFilePath(sourceId, "normalizeFile")
#     x <- readRDS(normalizeFilePath)
#     setkey(x$colData, "cell_id")
#     x
# }

# loadSampleCommon <- function(cacheKey, keyObject, key, cacheObject, sourceId, ...){
#     dprint("loadSampleCommon")
#     x <- loadSampleRds(sourceId)

#     x
#     # x[c("rowRanges", "colData")]
# }
# loadSampleWorking <- function(cacheKey, keyObject, key, cacheObject, sourceId, chrom){
#     dprint("loadSampleWorking")
#     x <- loadSampleRds(sourceId)
#     chrom_ <- chrom
#     isAll <- chrom == "all"
#     if(isAll){
#         goodBins <- x$rowRanges[, bin_n]
#     } else {
#         goodBins <- x$rowRanges[chrom == chrom_, chrom_bin_n]
#     }
#     nGoodBins <- length(goodBins)
#     allBins   <- 1:max(goodBins)
#     badBins   <- allBins[!(allBins %in% goodBins)]
#     nBadBins  <- length(badBins)    
#     nCells <- nrow(x$colData)
#     nullCell <- rep(NA,  nGoodBins) # we DO save all cells to allow dynamic cell rejection
#     nullMatrix <- matrix(NA, nrow = nBadBins, ncol = nCells)
#     order <- order(c(goodBins, badBins))
#     sd <- sapply(x$cells, function(cell){
#         i <- cell$hmm > 0
#         sd((cell$cn[i] - cell$hmm[i]) / cell$hmm[i], na.rm = TRUE) # NA if NULL, etc.
#     })
#     getZ <- function(val) lapply(-1:1, function(increment){
#         sapply(x$colData$cell_id, function(cell_id){
#             cell <- x$cells[[cell_id]]
#             if(is.null(cell[[val]])) return(nullCell)            
#             w <- x$windows[[paste("w", cell$window_size, sep = "_")]]
#             bins <- w[, if(isAll) TRUE else chrom == chrom_]
#             x <- if(val == "hmm"){
#                 ifelse( # hmm as RGB intensities
#                     cell[[val]][bins] == cell$modal_CN, 
#                     0.5,
#                     ifelse(
#                         sign(cell[[val]][bins] - cell$modal_CN) == sign(increment),
#                         1,
#                         0
#                     )
#                 )
#             } else { # cn as CNC z-scores
#                 cn <- cell$modal_CN + increment
#                 (cell[[val]][bins] - cn) / cn / sd[cell_id]                
#             }
#             uncollapseVector(x, cell$window_size, nGoodBins)    
#         })
#     })
#     padZ <- function(z) {
#         if(nBadBins > 0) lapply(z, function(m) rbind(m, nullMatrix)[order, ]) else z
#     }

#     valueTypes <- c("cn", "hmm")
#     windowMedians <- lapply(valueTypes, function(val){
#         d <- sapply(x$colData$cell_id, function(cell_id){
#             cell <- x$cells[[cell_id]]
#             if(is.null(cell[[val]])) return(nullCell)            
#             w <- x$windows[[paste("w", cell$window_size, sep = "_")]]
#             bins <- w[, if(isAll) TRUE else chrom == chrom_]
#             uncollapseVector(cell[[val]][bins], cell$window_size, nGoodBins)   
#         })
#         d <- apply(d, 1, median, na.rm = TRUE)
#         lapply(x$windows, function(w) {
#             n <- nGoodBins / w[chrom == chrom_, .N]
#             collapseVector(d, n) / n
#         })
#     })
#     names(windowMedians) <- valueTypes

#     getValues <- function(val){
#         x <- lapply(x$colData$cell_id, function(cell_id){
#             cell <- x$cells[[cell_id]]
#             if(is.null(cell[[val]])) return(NA)
#             wx <- paste("w", cell$window_size, sep = "_")
#             w <- x$windows[[wx]]
#             bins <- w[, if(isAll) TRUE else chrom == chrom_]
#             cell[[val]][bins] # - windowMedians[[val]][[wx]]
#         })
#         names(x) <- x$colData$cell_id
#         x
#     }
#     list(
#         cell_id     = x$colData$cell_id,
#         modal_CN    = x$colData$modal_CN,
#         rejected    = x$colData$rejected,
#         window_size = x$colData$window_size,
#         windows = lapply(x$windows, function(w) w[chrom == chrom_, .(start, end)]),
#         windowMedians = windowMedians,
#         cn  = getValues("cn"),
#         hmm = getValues("hmm"),
#         z = list(
#             cn  = padZ(getZ("cn")),
#             hmm = padZ(getZ("hmm"))
#         )
#     )
# }
