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
expandBatchValue <- function(x, windowPower, minPower, expandedWindows, force = FALSE){
    if(windowPower > minPower || force) {
        expandBy <- 2 ** (windowPower - minPower) # only expand to minPower, not necessarily individual bins
        if(force) expandBy <- 1 / expandBy
        x <- windows[[windowPower + 1]][, {
            chrom_ <- chrom
            nExpandedWindows <- expandedWindows[chrom == chrom_, .N]
            uncollapseVector(x[.I], expandBy, nExpandedWindows)
        }, by = chrom][[2]]
    } else if(windowPower < minPower){
        collapseBatchValue(x, windowPower, minPower, expandedWindows, force = TRUE) 
    }
    x
}
collapseBatchValue <- function(x, windowPower, minPower, expandedWindows, force = FALSE){
    if(windowPower > minPower || force) {
        expandBy <- 2 ** (windowPower - minPower)
        if(force) expandBy <- 1 / expandBy
        x <- windows[[minPower + 1]][, {
            collapseVector(x[.I], expandBy) / expandBy
        }, by = chrom][[2]]
    } else if(windowPower < minPower){
        expandBatchValue(x, windowPower, minPower, expandedWindows, force = TRUE)
    }
    x
}
getBatchWindows <- function(cell){
    repKey <- if(cell$cellIsReplicating) "composite" else "sequential"
    cell$windows[[shapeKey]][[repKey]] 
}
getExpandedHMM <- function(cells, minPower, expandedWindows, maskExtremes = TRUE){
    byCell <- simplify2array(mclapply(cells, function(cell){
        cellWindows <- getBatchWindows(cell)
        expandBatchValue(cellWindows$HMM, cell$windowPower, minPower, expandedWindows)
    }, mc.cores = env$N_CPU))
    windowMedians <- round(apply(byCell, 1, median, na.rm = TRUE), 0)
    if(maskExtremes) windowMedians <- ifelse(windowMedians > 0 & windowMedians < 4, windowMedians, NA)
    list(
        byCell = byCell,
        median = windowMedians
    )
}
#=====================================================================================

#=====================================================================================
# main target function that launches the inter-sample batch normalization process
#-------------------------------------------------------------------------------------
normalizeBatch <- function(sampleName, cells){

    # get and check for sufficient working cells
    cell_ids <- sapply(cells, function(cell) cell$cell_id)
    cell_ids_working <- colData[cell_id %in% cell_ids & !bad & keep, cell_id]
    nCells <- length(cell_ids)
    nWorkingCells <- length(cell_ids_working)
    if(nWorkingCells < 5){
        message(paste("  sample", sampleName, ": only", nWorkingCells, "of", nCells, "cells are usable, skipping batch correction"))
        return(cells)
    }
    message(paste("  sample", sampleName, ":", nWorkingCells, "of", nCells, "cells used for batch fit"))
    workingCells <- cells[cell_ids_working]

    # determine the common, expanded window basis for performing cell comparisons
    allCellWindowPowers <- sapply(cells, function(cell) cell$windowPower)
    workingCellWindowPowers <- sapply(workingCells, function(cell) cell$windowPower)
    allWindowPowers <- sort(unique(allCellWindowPowers))
    workingWindowPowers <- sort(unique(workingCellWindowPowers))
    minWorkingWindowPower <- min(workingWindowPowers)
    expandedWindows <- windows[[minWorkingWindowPower + 1]] # the highest resolution windows used by the sample's cells
    
    # determine the median HMM call per expanded window; take as the sample's reference, baseline CN
    expandedHMM <- getExpandedHMM(workingCells, minWorkingWindowPower, expandedWindows)
    
    # determine the CN standard deviation of each cell normalized to CN == 1
    # do this for all cells, not just those being used for batch correction
    cncsd_norm <- mclapply(cells, function(cell){
        if(cell$badCell) return(NA)
        cellWindows <- getBatchWindows(cell)
        I <- cellWindows$HMM == cell$ploidy # restrict SD calculation to the most abundant, highest reliability windows
        sd((cellWindows$CN[I] - cell$ploidy), na.rm = TRUE) / sqrt(cell$ploidy)   
    }, mc.cores = env$N_CPU)
    names(cncsd_norm) <- cell_ids

    # establish CN Z scores for every expanded window for each working cell
    zScores <- simplify2array(mclapply(workingCells, function(cell){
        cellWindows <- getBatchWindows(cell)
        CN <- expandBatchValue(cellWindows$CN, cell$windowPower, minWorkingWindowPower, expandedWindows)
        ifelse( # restrict Z score calculation to the sample HMM reference, again, the most reliable
            expandedHMM$byCell[, cell$cell_id] == expandedHMM$median,
            (CN - expandedHMM$median) / sqrt(expandedHMM$median) / cncsd_norm[[cell$cell_id]], 
            NA
        )
    }, mc.cores = env$N_CPU))
    
    # establish a weighting scheme when aggregating windows from cells with different windowPowers
    windowPowerPairs <- as.data.table(expand.grid(allWindowPowers, allWindowPowers))
    setnames(windowPowerPairs, c("wp1", "wp2"))
    windowPowerPairs$key <- apply(windowPowerPairs, 1, paste, collapse = ":")         
    windowPowerPairs[, weight := 2 ** pmin(wp1, wp2) / 2 ** pmax(wp1, wp2)]        
    setkey(windowPowerPairs, "key")

    # calculate the weighted median Z score for expanded windows for all cell windowPowers in use
    medianZScores <- mclapply(allWindowPowers, function(wp1){
        keys <- paste(wp1, workingCellWindowPowers, sep = ":") # cannot embed directly in windowPowerPairs[]
        weights <- windowPowerPairs[keys, weight]
        apply(zScores, 1, function(x) weighted.median(x, weights))
    }, mc.cores = env$N_CPU)
    names(medianZScores) <- as.character(allWindowPowers)

    # rescale median Z scores back to copy number units for each cell
    maxBatchShape <- 3
    mclapply(cells, function(cell){
        if(cell$badCell) return(cell)
        cellWindows <- getBatchWindows(cell)
        shape <- (1 + medianZScores[[as.character(cell$windowPower)]] * cncsd_norm[[cell$cell_id]]) / 1
        shape <- collapseBatchValue(shape, cell$windowPower, minWorkingWindowPower, expandedWindows)     
        cell$windows$batched <- list(
            shape = pmin(maxBatchShape, pmax(1 / maxBatchShape, shape))
        )
        cell
    }, mc.cores = env$N_CPU)
}
#=====================================================================================

#=====================================================================================
# regardless of whether cells were batch corrected, calculate the median HMM and sex
#-------------------------------------------------------------------------------------
getSampleProfile <- function(cells){

    # get and check for sufficient working cells
    cell_ids <- sapply(cells, function(cell) cell$cell_id)
    cell_ids_working <- colData[cell_id %in% cell_ids & !bad & keep, cell_id]
    if(length(cell_ids_working) == 0) return(NA)    
    workingCells <- cells[cell_ids_working]

    # determine the median HMM call per expanded window; take as the sample's reference, baseline CN
    expandedHMM <- getExpandedHMM(workingCells, 0, windows[[1]], maskExtremes = FALSE)

    # provide the cell reference HMM CN for all windowPowers
    mclapply(windowPowers, function(windowPower){
        round(collapseBatchValue(expandedHMM$median, windowPower, 0, windows[[1]]), 0)
    }, mc.cores = env$N_CPU)
}
collateChromosomes <- function(windows, cnProfile){

    # find the median HMM per chromosome, take it as the chromosome's copy number
    CN <- windows[, 
        .(
            observed = round(median(cnProfile[.I], na.rm = TRUE), 0),
            expected = env$PLOIDY
        ), 
        by = "chrom"
    ]

    # is this a typical mammalian genome with XY sex chromosomes?
    setkey(CN, chrom)
    sexChroms <- c("chrX", "chrY") 
    if(CN[, sum(chrom %in% sexChroms) != 2]) return(list(CN = CN, sex = as.character(NA)))

    # if so, determine the apparent sex as XX, XY, XO, XYY, etc.
    CN_X <- CN["chrX"]$observed
    CN_Y <- CN["chrY"]$observed
    X <- paste0(rep("X", CN_X), collapse = "")
    Y <- paste0(rep("Y", CN_Y), collapse = "")
    if(env$PLOIDY == 2){
        if(CN_X == 0) X <- "O" 
        if(CN_X < 2 && CN_Y == 0) Y <- "O"
        if(CN_X < 2){ # expect male
            expectedX <- 1
            expectedY <- 1
        } else { # expect female
            expectedX <- 2
            expectedY <- 0
        }
    } else {
        if(CN_X > 0 && CN_Y == 0){ # expect X sperm
            expectedX <- 1
            expectedY <- 0
        } else if(CN_X == 0 && CN_Y > 0){ # expect Y sperm
            expectedX <- 0
            expectedY <- 1
        } else { # for XY sperm (or other haploid) impossible to set a single expectation
            expectedX <- 0
            expectedY <- 0
        }
    }
    CN[chrom == "chrX", expected := expectedX] 
    CN[chrom == "chrY", expected := expectedY] 
    return(list(CN = CN, sex = paste0(X, Y)))
}
#=====================================================================================
