# determine the least number of reads per bin for a cell's data to support a robust HMM
# places ~96% of bins within modal_CN +/- 0.5
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

# demarcation lines for plots
NR_map_hLines <- function(cellData, cellWindowCounts){

    modal_CN <- cellData()$colData$modal_CN     

    ER_map_modal_CN <- cellWindowCounts()$ER_map_modal_CN
    readsPerAllele <- ER_map_modal_CN / modal_CN
    orange <- CONSTANTS$plotlyColors$orange
    red    <- CONSTANTS$plotlyColors$red
    structure(
        ER_map_modal_CN + -2:2 * readsPerAllele, 
        color = c(red, red, orange, red, red)
    )
}
