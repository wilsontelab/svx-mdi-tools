#----------------------------------------------------------------------
# handle binned coverage loading and filtering
# these are generic functions that apply to potentially any SVX app
# expects:
#   app-specific loadFn() that returns
# list(
#     maxPower = 1,
#     binDensities = sapply(0:maxPower, function(power) binDensity0 / aggFactor**power), # aggFactor typically 10
#     binSizes     = sapply(0:maxPower, function(power) binSize0    * aggFactor**power),
#     medianCoverage = medianCoverage,
#     x0 = x0[, .(chrom,start,genomeStart,coverage)], # power 0 (i.e., as is)
#     x1 = x1[, .(chrom,start,genomeStart,coverage)]  # power 1 (i.e., 10-fold pre-merging of bins), etc.
# )
#----------------------------------------------------------------------

# pre-calculate aggregated bins to speed up window loads at low resolution
svx_rebinCoverage <- function(sourceId, sample, x0, aggFactor = 10){
    startSpinner(session, message = "rebinning coverage .")
    binSize <- x0[1:2, diff(start)]
    binDensity <- 1 / binSize
    # x0[, ":="(
    #     index1 = floor(genomeStart / binSize / aggFactor) # for the next round aggregation
    # )]

    # startSpinner(session, message = "rebinning coverage ..")
    # x1 <- x0[, .(
    #     start = start[1], 
    #     genomeStart = genomeStart[1],            
    #     coverage = mean(coverage, na.rm = TRUE),
    #     gc = mean(gc),
    #     excluded = mean(excluded),
    #     index2 = floor(index1 / aggFactor)
    # ), by = .(chrom, index1)]

    # startSpinner(session, message = "rebinning coverage ...")
    # x2 <- x1[, .(
    #     start = start[1],
    #     genomeStart = genomeStart[1],            
    #     coverage = mean(coverage, na.rm = TRUE),
    #     gc = mean(gc),
    #     excluded = mean(excluded)
    # ), by = .(chrom, index2)]

    startSpinner(session, message = "fitting coverage")
    # medianCoverage <- x2[, {
    medianCoverage <- x0[, {
        median <- median(coverage, na.rm = TRUE)
        peakValue(coverage[coverage > median / 5 & coverage < median * 5])
    }]

    startSpinner(session, message = "aggregating coverage")
    outCols <- c("chrom","start","genomeStart","coverage","gc","excluded")    
    # maxPower <- 2
    maxPower <- 0
    list(
        sourceId = sourceId,
        sample = sample,
        maxPower = maxPower,
        binDensities = sapply(0:maxPower, function(power) binDensity / aggFactor**power),
        binSizes     = sapply(0:maxPower, function(power) binSize    * aggFactor**power),
        medianCoverage = medianCoverage,
        x0 = x0[, .SD, .SDcols = outCols]
        # ,
        # x1 = x1[, .SD, .SDcols = outCols],
        # x2 = x2[, .SD, .SDcols = outCols]
    )
}

# filter coverage profile by browser coordinates
svx_filterCoverageByRange <- function(sourceId, sample, coord, maxBins, loadFn){
    coverage <- loadFn(sourceId, sample)  
    startSpinner(session, message = "filtering coverage")
    # maxBinDensity <- maxBins / coord$width
    # passingPowers <- which(coverage$binDensities <= maxBinDensity) - 1L
    # power <- if(length(passingPowers) == 0) coverage$maxPower
    #     else if(passingPowers[1] == 0) 0 
    #     else passingPowers[1] - 1 # thus, if possible, return the first density that did NOT pass, to allow partial aggregation below
    power <- 0
    x <- coverage[[paste0("x", power)]] 
    isWholeGenome <- coord$chromosome == "all"
    isProperChrom <- isProperChromosome(coord$chromosome)
    isGenome <- isWholeGenome || !isProperChrom
    x <- if(isGenome) x[data.table::between(as.numeric(genomeStart), as.numeric(coord$start), as.numeric(coord$end))]
                 else x[data.table::between(start,                   as.numeric(coord$start), as.numeric(coord$end)) & chrom == coord$chromosome]
    list(
        sourceId = coverage$sourceId,
        sample   = coverage$sample,
        binSize  = coverage$binSizes[power + 1],
        medianCoverage = coverage$medianCoverage,
        bins = x[, .(
            chrom,
            gc,
            excluded,
            strand = ".",
            x =if(isGenome) genomeStart else start,
            y = coverage
        )]
    )
}

# support read depth and copy number plots
svx_setCoverageValue <- function(coverage, plotAs, medianPloidy){
    coverage$bins[, .(
        strand = ".",
        x = x,
        y = switch(
            plotAs,
            "Read Depth" = y,
            {
                normalized <- if(is.null(app$normalizeGC)) NULL 
                              else app$normalizeGC$getBinNormalizedCN(coverage$sourceId, coverage$sample, "uncollapsed", gc, y)
                if(is.null(normalized)) y / coverage$medianCoverage * medianPloidy else normalized
            }
        ),
        z = excluded
    )]
}

# mask low quality bins from coverage plots
svx_maskLowQualityBins <- function(bins){
    bins[z <= 0.1, .SD, .SDcols = c("strand", "x", "y")]
}

# collect copy number data from normalizeGC, if in use
svx_getCnvJxnNormalizedCN <- function(isMultiSample, targetId){
    if(!isMultiSample || is.null(app$normalizeGC)) return(list(key = NA, value = NA))
    app$normalizeGC$getCnvJxnsNormalizedCN(targetId)[c("key","value")]
}
svx_getCnvJxnNormalizedCN_singleJunction <- function(targetId, svId){
    if(is.null(app$normalizeGC)) return(NA)
    app$normalizeGC$getCnvJxnsNormalizedCN(targetId)$value$dt[SV_ID == svId, maxCNC]
}

svx_getHmmCnvs <- function(targetId){
    if(is.null(app$normalizeGC)) return(list(key = NA, value = NA))
    app$normalizeGC$getNormalizedHmmCnvs(targetId)[c("key","value")]
}