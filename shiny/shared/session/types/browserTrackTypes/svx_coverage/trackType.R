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
        I <- sapply(names(selectedSources), function(x) sampleName %in% selectedSources[[x]]$Sample_ID)
        sourceId <- names(selectedSources)[I] 
        x <- svx_filterCoverageByRange(sourceId, sampleName, coord, maxBins, loadFn) %>%
        svx_setCoverageValue(reference, coord, Plot_Type, Median_Ploidy) %>%
        aggregateTabixBins(track, coord, plotBinSize) %>%
        svx_maskLowQualityBins()
    }
    cnvHighlightsFn <- function(track, reference, coord, sampleName, sample){
        I <- sapply(names(selectedSources), function(x) sampleName %in% selectedSources[[x]]$Sample_ID)
        x <- svx_getHmmCnvs(names(selectedSources)[I])
        if(!isTruthy(x)) return(NULL)
        x$value[[sampleName]][
            chrom == coord$chromosome &
            start <= coord$end & 
            coord$start <= end &
            !is.na(JXN_TYPE),
            .(
                x1 = start,
                x2 = end,
                color = ifelse(JXN_TYPE == "D", rgb(1,0.2,0.2,0.075), rgb(0.2,0.2,1,0.075))
            )
        ]
    }

    # spans <- sapply(names(selectedSources), svx_getHmmCnvs, simplify = FALSE, USE.NAMES = TRUE)
    # dstr(spans)
# List of 1
#  $ e59b0384bfa375dc7846d5c384d62450:List of 2
#   ..$ key  : chr "7bebd41b95d0709bcd12882c378e014f"
#   ..$ value:List of 2
#   .. ..$ HCT116 :Classes ‘data.table’ and 'data.frame': 410 obs. of  13 variable
# s:

#   .. .. ..$ JXN_TYPE  : chr [1:410] "L" "L" "L" "L" ...


    # build the binned_XY_track
    isDifference <- nSamples > 1 && getTrackSetting(track, "Data", "Aggregate", "none") == "difference"
    build.binned_XY_track(
        track, reference, coord, layout, 
        dataFn, 
        stranded = FALSE, 
        allowNeg = isDifference, 
        ylab = if(isDifference) paste(Plot_Type, "Change") else Plot_Type,
        center = TRUE, binSize = plotBinSize
        # ,
        # highlightsFn = cnvHighlightsFn
    )
}
