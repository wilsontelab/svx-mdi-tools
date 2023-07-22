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
build.svx_coverageTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)

    Max_Points          <- getBrowserTrackSetting(track, "Coverage", "Max_Points", 1000)
    Max_Copy_Number     <- getBrowserTrackSetting(track, "Coverage", "Max_Copy_Number", 6) 
    Max_Coverage        <- getBrowserTrackSetting(track, "Coverage", "Max_Coverage", 100) 
    Median_Ploidy       <- getBrowserTrackSetting(track, "Coverage", "Median_Ploidy", 2)
    Plot_As             <- getBrowserTrackSetting(track, "Coverage", "Plot_As", "Read Depth")  
    Combine_Samples     <- getBrowserTrackSetting(track, "Coverage", "Combine_Samples", FALSE) 





    sampleItems <- track$settings$items()
    selectedSources <- getSourcesFromTrackSamples(sampleItems)
    selectedSamples <- list()
    sampleI <- 1



    
    coverage <- do.call(rbind, lapply(names(selectedSources), function(sourceId){
        do.call(rbind, lapply(selectedSources[[sourceId]]$Sample_ID, function(sample_){
            selectedSamples[[sample_]] <<- sampleI
            sampleI <<- sampleI + 1
            filterCoverageByRange(sourceId, sample_, coord, Max_Points) %>% 
            aggregateCoverageBins(Max_Points) %>% 
            setCoverageXY(Plot_As, Median_Ploidy, coord, sample_)
        }))
    }))

    if(Plot_As == "Copy Number Change"){
        req(length(selectedSamples) > 1)
        refSample <- sampleItems[[1]]$Sample_ID
        refY <- coverage[sample == refSample, y]
        coverage <- coverage[sample != refSample, .(
            x = x,
            y = y - refY
        ), by = .(sample)][, .SD, .SDcols = c("x","y","sample")]    
        selectedSamples <- selectedSamples[2:length(selectedSamples)]      
    }

    if(Combine_Samples){
        selectedSamples <- list(combined = 1)
        coverage <- coverage[, .(
            y = if(Plot_As == "Read Depth") sum(y) else mean(y),
            sample = "combined"
        ), by = .(x)]
    }

    ylim <- switch(
        Plot_As,
        'Read Depth' = {
            if(Max_Coverage == 0) Max_Coverage <- max(coverage$y) * 1.05
            c(0, Max_Coverage)
        },
        'Copy Number' = {
            if(Max_Copy_Number == 0) Max_Copy_Number <- max(coverage$y) * 1.05
            c(0, Max_Copy_Number)
        },
        'Copy Number Change' = {
            if(Max_Copy_Number == 0) Max_Copy_Number <- max(abs(coverage$y)) * 1.05
            c(-Max_Copy_Number, Max_Copy_Number)
        }
    ) 

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches    
    image <- mdiTrackImage(layout, height, message = getBrowserTrackSetting(track, "Track_Options", "Track_Name", "svx_coverage"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))

        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = Plot_As, #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        abline(
            h = if(Plot_As == "Read Depth") seq(0, ylim[2], 10) else -Max_Copy_Number:Max_Copy_Number, 
            col = CONSTANTS$plotlyColors$grey
        )

        I <- if(length(selectedSamples) > 1) coverage[, sample(.N)] else TRUE
        colI <- coverage[I, unlist(selectedSamples[sample])]
        points(coverage[I], pch = 19, col = unlist(CONSTANTS$plotlyColors[colI]))

        colI <- unlist(selectedSamples)
        trackLegend(
            track, coord, ylim, 
            legend = names(selectedSamples), pch = 19, col =  unlist(CONSTANTS$plotlyColors[colI])
        )
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}
