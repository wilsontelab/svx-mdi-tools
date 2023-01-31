# determine the least number of reads per bin for a cell's data to support a robust HMM
# places ~96% of bins within modal_CN +/- 0.5
# sapply(1:4, getMinBinCount) => 16  64 144 256; thus, usually 64 for diploid cells
getMinBinCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
    (CONSTANTS$nSdHalfCn / (ploidyFactor - 1)) ** 2
}

# set the per-cell window size as the number of bins needed to obtain a mean raw count >=minBinCount
getMinWindowPower <- function(cellData){
    minBinCount <- getMinBinCount(cellData$colData[, modal_CN])
    NR_raw_b <- cellData$rawCounts
    NR_avg_b <- mean(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_avg_b) || NR_avg_b == 0) return(CONSTANTS$maxWindowPower)
    window_size <- ceiling(minBinCount / NR_avg_b) # in number of bins (not bp)
    if(is.na(window_size) || window_size > CONSTANTS$maxWindowSize) return(CONSTANTS$maxWindowPower)
    window_size <- CONSTANTS$windowSizes[min(which(CONSTANTS$windowSizes >= window_size))]
    log2(window_size)
}

# TODO: move this to pipeline as new approach to establishing cell parameters
# can then pre-fit GC bias to the default window size, which is generally reliable
# but before doing the fit, should attempt to find:
#   large segmental changes, including aneuploidy - HMM achieves this already
#   replicating cells, via a mixture model 
#   do this only on segments at modal CN, i.e., consistent with majority genome

# use bin-to-bin count deltas to determine the window size that
# puts most window counts between CN = ploidy +/- 1
# tentatively keep cells that can achieve this, reject cells that can't
getCellDefaults <- function(mappableWindows, cellData){
    getSdLagDiff <- function(windowPower){
        mappableWindows <- mappableWindows[[windowPower + 1]]
        windowSize <- as.integer(2 ** windowPower)    
        NR <- collapseVector(cellData$rawCounts, windowSize)[mappableWindows$i]
        NR_map <- NR / mappableWindows$windows$mappability
        NR_map <- NR_map[!is.na(NR_map) & NR_map > 0]
        lagDiff <- diff(NR_map)
        readPerAllele <- peakValue(NR_map) / cellData$colData$modal_CN
        sd(lagDiff) / readPerAllele * cellData$colData$modal_CN        
    }
    minWindowPower <- getMinWindowPower(cellData)
    windowPower <- minWindowPower
    sdLagDiff <- getSdLagDiff(windowPower)
    while(sdLagDiff > 1 && windowPower < CONSTANTS$maxWindowPower){
        windowPower <- windowPower + 1
        sdLagDiff <- getSdLagDiff(windowPower)
    }
    list(
        windowPower = windowPower,
        keep = sdLagDiff <= 1
    )
}

# demarcation lines for plots
NR_map_hLines <- function(modal_CN, cellWindowCounts){
    ER_map_modal_CN <- cellWindowCounts$ER_map_modal_CN
    readsPerAllele <- ER_map_modal_CN / modal_CN
    orange <- CONSTANTS$plotlyColors$orange
    red    <- CONSTANTS$plotlyColors$red
    structure(
        ER_map_modal_CN + -5:5 * readsPerAllele, 
        color = c(rep(red, 5), orange, rep(red, 5))
    )
}

# node handling
addNode <- function(session, plotId, x, y){
    plotId <- paste(plotId, "plotly", sep ="-")
    p <- plotlyProxy(plotId, session, deferUntilFlush = FALSE)   
    plotlyProxyInvoke(p, "extendTraces", list(
        x = list(list(x)),
        y = list(list(y))
    ), list(1)) # 1 is the index of the overplot trace 
}
clickPlotNode <- function(session, plot, plotId, nodes){
    d <- plot$clicked()
    matchingNode <- nodes[, x == d$x]
    if(!any(matchingNode)) {
        addNode(session, plotId, d$x, d$y)
        rbind(nodes, data.table(x = d$x, y = d$y)) 
    }
}
