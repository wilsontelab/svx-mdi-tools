# determine the least number of reads per bin for a cell's data to support a robust HMM
# places ~96% of bins within modal_CN +/- 0.5
N_SD_HALFCN <- 2
getMinBinCount <- function(modal_CN){
    ploidyFactor <- (modal_CN + 0.5) / modal_CN
    (N_SD_HALFCN / (ploidyFactor - 1)) ** 2
}

# set the per-cell window size as the number of bins needed to obtain a mean raw count >=minBinCount
MAX_WINDOW_SIZE <- 64
MAX_WINDOW_POWER <- log2(MAX_WINDOW_SIZE)
WINDOW_SIZES <- 2 ** (0:MAX_WINDOW_POWER)
getMinWindowSize <- function(sourceData, cell_id){
    minBinCount <- getMinBinCount(sourceData$colData[cell_id, modal_CN])
    NR_raw_b <- sourceData$raw_counts[[cell_id]]
    NR_avg_b <- mean(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_avg_b) || NR_avg_b == 0) return(MAX_WINDOW_SIZE)
    window_size <- ceiling(minBinCount / NR_avg_b)
    if(is.na(window_size) || window_size > MAX_WINDOW_SIZE) return(MAX_WINDOW_SIZE)
    WINDOW_SIZES[min(which(WINDOW_SIZES >= window_size))]
}

# demarcation lines for plots
getChromLines <- function(cellWindows){
    L <- cellWindows()$lastChromBins
    structure(c(0, L) + 0.5, color = rep('grey', length(L) + 1))    
}
getPeakNR <- function(cellWindows){
    v <- peakValue(cellWindows()$NR_map)
    structure(v, color = 'red')
}
