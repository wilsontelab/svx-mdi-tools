#----------------------------------------------------------------------
# svx_coverage trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_svx_coverageTrack <- function(trackId) {
    list( 
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        expand = NULL,
        NULL
    )
}

# build method for the S3 class; REQUIRED
build.svx_coverageTrack <- function(track, reference, coord, layout, loadFn){
    req(coord, coord$chromosome)

    # get operating parameters
    Plot_Type       <- getTrackSetting(track, "Coverage", "Plot_Type", "Read Depth")  
    Median_Ploidy   <- getTrackSetting(track, "Coverage", "Median_Ploidy", 2)
    plotBinSize     <- as.integer(getTrackSetting(track, "X_Axis", "Bin_Size", ""))
    if(!isTruthy(plotBinSize)) {
        maxBins <- getTrackSetting(track, "X_Axis", "Max_Bins", 250)
        plotBinSize <- floor(coord$width / maxBins)
    } else {
        maxBins <- floor(coord$width / plotBinSize)
    }

    # parse our binned coverage data source
    nSamples <- length(track$settings$items())
    selectedSources <- getSourcesFromTrackSamples(track$settings$items())
    dataFn <- function(track, reference, coord, sampleName, sample){
        sourceId <- c(sapply(names(selectedSources), function(x) if(sampleName %in% selectedSources[[x]]$Sample_ID) x else NULL))
        x <- svx_filterCoverageByRange(sourceId, sampleName, coord, maxBins, loadFn)
        aggregateTabixBins(x$bins, track, coord, plotBinSize) %>%
        svx_setCoverageValue(Plot_Type, x$medianCoverage, Median_Ploidy)
    }

    # build the binned_XY_track
    isDifference <- nSamples > 1 && getTrackSetting(track, "Data", "Aggregate", "none") == "difference"
    build.binned_XY_track(
        track, reference, coord, layout, 
        dataFn, 
        stranded = FALSE, 
        allowNeg = isDifference, 
        ylab = if(isDifference) paste(Plot_Type, "Change") else Plot_Type,
        center = TRUE, binSize = plotBinSize
    )
}
