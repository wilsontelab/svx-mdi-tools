#=====================================================================================
# use a batch of cells to normalize away consistent bin deviations from expectations
# note! this includes CNVs commons to the majority of cells!
#=====================================================================================

#=====================================================================================
# plotting functions
#-------------------------------------------------------------------------------------
saveBatchPlot <- function(name, fn, width = 2, height = 2){ # save a plot for future assembly and app
    if(!env$VERBOSE_PLOTS) return(NULL)
    pngFile <- file.path(plotsDir, paste0(name, ".png"))
    png(pngFile, width = width, height = height, units = "in", pointsize = 7, res = 300, type = "cairo")
    par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
    fn()
    dev.off()
}
#=====================================================================================

#=====================================================================================
# support functions
#-------------------------------------------------------------------------------------
weighted.median <- function(x, w) { # adapted from spatstat
    I <- !(is.na(x) | is.na(w))
    x <- x[I]
    w <- w[I]
    if(length(x) == 0) return(NA)
    if(length(x) == 1) return(x)
    xOrder <- order(x)
    x <- x[xOrder]
    w <- w[xOrder]
    Fx <- cumsum(w) / sum(w)
    if(anyDuplicated(x)) {
        dup <- rev(duplicated(rev(x)))
        x  <- x[!dup]
        Fx <- Fx[!dup]
    }
    approx(Fx, x, xout = 0.5, ties = "ordered", 
           rule = 2, method = "constant", f = 1/2)$y
}
# NR_wmB_buffer <- list()
# getNR_wmB <- function(cell) with(cell, { # apply the batch shape correction to NR_wm
#     if(!is.list(batchShapes) || env$SHAPE_CORRECTION != "batch") return(windows$NR_wm)
#     if(is.null(NR_wmB_buffer[[cell_id]])) {
#         NR_wmB <- windows$NR_wm / batchShapes[[as.character(windowPower)]]       
#         N_pre  <- sum(windows$NR_wm, na.rm = TRUE)
#         N_post <- sum(NR_wmB, na.rm = TRUE)
#         NR_wmB_buffer[[cell_id]] <<- NR_wmB * N_pre / N_post
#     }
#     NR_wmB_buffer[[cell_id]]
# })
#=====================================================================================

#=====================================================================================
# main target function that launches the inter-sample batch normalization process
#-------------------------------------------------------------------------------------
normalizeBatch <- function(cells){

    # determine the common, expanded window basis for performing cell comparisons
    cellWindowPowers <- sapply(cells, function(cell) cell$windowPower)
    workingWindowPowers <- sort(unique(cellWindowPowers))
    minCellWindowPower <- min(workingWindowPowers)
    expandedWindows <- windows[[minCellWindowPower + 1]] # the highest resolution windows used by the sample's cells

    # determine the median HMM CN call per expanded window; take as the sample standard
    expandValue <- function(x, cell){
        if(cell$windowPower > minCellWindowPower) {
            windows <- windows[[cell$windowPower + 1]]
            expandBy <- 2 ** (cell$windowPower - minCellWindowPower) # only expand to minCellWindowPower, not necessarily individual bins
            x <- windows[, {
                chrom_ <- chrom
                nExpandedWindows <- expandedWindows[chrom == chrom_, .N]
                uncollapseVector(x[.I], expandBy, nExpandedWindows)
            }, by = chrom][[2]]
        }
        x
    }
    collapseValue <- function(){

    }
    expandedHMM <- sapply(cells, function(cell){
        cellWindows <- cell$windows[[cell$windowPower + 1]][[shapeKey]]$sequential
        expandValue(cellWindows$HMM, cell)
    })
    medianExpandedHMM <- round(apply(expandedHMM, 1, median, na.rm = TRUE), 0)

    # determine the CN standard deviation of each cell
    cell_ids <- sapply(cells, function(cell) cell$cell_id)
    cncsd_norm <- lapply(cells, function(cell){
        cellWindows <- cell$windows[[cell$windowPower + 1]][[shapeKey]]$sequential 
        I <- cellWindows$HMM == cell$ploidy # restrict SD calculation to the most abundant, highest reliability windows
        sd((cellWindows$CN[I] - cell$ploidy), na.rm = TRUE) / sqrt(cell$ploidy)   
    })
    names(cncsd_norm) <- cell_ids

    # establish CN Z scores for every expanded window for each cell
    zScores <- sapply(cells, function(cell){
        cellWindows <- cell$windows[[cell$windowPower + 1]][[shapeKey]]$sequential
        CN  <- expandValue(cellWindows$CN,  cell)
        ifelse( # restrict Z score calculation to the sample HMM standard, again, the most reliable
            expandedHMM[, cell$cell_id] == medianExpandedHMM, 
            (CN - medianExpandedHMM) / sqrt(medianExpandedHMM) / cncsd_norm[[cell$cell_id]], 
            NA
        )
    })
    
    # establish a weighting scheme when aggregating windows from cells with different windowPowers
    windowPowerPairs <- as.data.table(expand.grid(workingWindowPowers, workingWindowPowers))
    setnames(windowPowerPairs, c("wp1", "wp2"))
    windowPowerPairs$key <- apply(windowPowerPairs, 1, paste, collapse = ":")         
    windowPowerPairs[, weight := 2 ** pmin(wp1, wp2) / 2 ** pmax(wp1, wp2)]        
    setkey(windowPowerPairs, "key")

    # calculate the weighted median Z score for expanded windows for all cell windowPowers in use
    medianZScores <- lapply(workingWindowPowers, function(wp1){
        keys <- paste(wp1, cellWindowPowers, sep = ":") # cannot embed directly in windowPowerPairs[]
        weights <- windowPowerPairs[keys, weight]
        # expandBy <- 2 ** (wp1 - minCellWindowPower)
        medianZScore <- apply(zScores, 1, function(x){
            if(all(is.na(x))) return(NA)
            weighted.median(x, weights)
            # if(x < 0.1) 1 else x # don't reshape bins with an effectively zero median
        })
        # if(expandBy > 1) medianZScore <- windows[[minCellWindowPower + 1]][, {
        #     collapseVector(medianZScore[.I], expandBy) / expandBy
        # }, by = chrom][[2]]
        medianZScore
    })
    names(medianZScores) <- as.character(workingWindowPowers)

    saveBatchPlot("medianZScores", function(){
        plot(density(medianZScores[[1]], na.rm = TRUE))
        # plot(1:length(medianZScore), log2(medianZScore), pch = 19, cex = 0.3, ylim = c(-3, 3),
        #     xaxt = "n", xaxs = 'i', xlab = "Genome Window", ylab = "log2(Shape)", col = defaultPointColor)
        abline(v = -3:3, col = "grey")
    }, width = 6)

    # rescale median Z scores back to copy number units for each cell
    batchShapes <- sapply(cells, function(cell){
        shape <- (1 + medianZScores[[as.character(cell$windowPower)]] * cncsd_norm[[cell$cell_id]]) / 1
        expandBy <- 2 ** (cell$windowPower - minCellWindowPower)
        print(cell$cell_id)
        print(cell$windowPower)
        print(expandBy)
        str(shape) 
        str(windows[[cell$windowPower + 1]])
        if(expandBy > 1) shape <- windows[[minCellWindowPower + 1]][, {
            collapseVector(shape[.I], expandBy) / expandBy
        }, by = chrom][[2]]
        str(shape)            
        pmin(10, pmax(0.1, shape))
    })

    saveBatchPlot("batchShapes", function(){
        plot(density(batchShapes[[1]], na.rm = TRUE))
        abline(v = 0:2, col = "grey")
    }, width = 6)

    # str(expandedHMM)
    # str(medianExpandedHMM)
    # str(cncsd_norm)
    # str(zScores)
    # str(medianZScores)
    # print(range(medianZScores, na.rm = TRUE))
    str(batchShapes)
    print(range(unlist(batchShapes), na.rm = TRUE))
    # str(windows)

    cellI <- 2
    saveBatchPlot("before", function(){
        cell <- cells[[cellI]]
        cellWindows <- cell$windows[[cell$windowPower + 1]][[shapeKey]]$sequential
        plot(1:length(cellWindows$CN), cellWindows$CN, pch = 19, cex = 0.3, col = defaultPointColor, ylim = c(0, 5))
        abline(h = 0:5, col = "grey")
    }, width = 6)
    saveBatchPlot("after", function(){
        cell <- cells[[cellI]]
        cellWindows <- cell$windows[[cell$windowPower + 1]][[shapeKey]]$sequential        
        str(cellWindows$CN)
        str(batchShapes[[cell$cell_id]])

        plot(1:length(cellWindows$CN), cellWindows$CN / batchShapes[[cell$cell_id]], pch = 19, cex = 0.3, col = defaultPointColor, ylim = c(0, 5))
        abline(h = 0:5, col = "grey")
    }, width = 6)

    batchShapes
}
#=====================================================================================
