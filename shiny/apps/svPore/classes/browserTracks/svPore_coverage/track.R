#----------------------------------------------------------------------
# svPore_coverage trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_svPore_coverageTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        expand = NULL,
        NULL
    )
}

# build method for the S3 class; REQUIRED
build.svPore_coverageTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)

    Max_Points          <- getBrowserTrackSetting(track, "Coverage", "Max_Points", 1000)
    Max_Copy_Number     <- getBrowserTrackSetting(track, "Coverage", "Max_Copy_Number", 6) 
    Max_Coverage        <- getBrowserTrackSetting(track, "Coverage", "Max_Coverage", 100) 
    Median_Ploidy       <- getBrowserTrackSetting(track, "Coverage", "Median_Ploidy", 2)
    Plot_As             <- getBrowserTrackSetting(track, "Coverage", "Plot_As", "Read Depth")  
    Combine_Samples     <- getBrowserTrackSetting(track, "Coverage", "Combine_Samples", FALSE) 

    sampleItems <- track$settings$items()
    sourcesToPlot <- getSvPoreSampleSources(sampleItems)
    samplesToPlot <- list()
    sampleI <- 1
    coverage <- do.call(rbind, lapply(names(sourcesToPlot), function(sourceId){
        do.call(rbind, lapply(sourcesToPlot[[sourceId]]$Sample_ID, function(sample_){
            samplesToPlot[[sample_]] <<- sampleI
            sampleI <<- sampleI + 1
            filterCoverageByRange(sourceId, sample_, coord, Max_Points) %>% 
            aggregateCoverageBins(Max_Points) %>% 
            setCoverageXY(Plot_As, Median_Ploidy, coord, sample_)
        }))
    }))

    if(Plot_As == "Copy Number Change"){
        req(length(samplesToPlot) > 1)
        refSample <- sampleItems[[1]]$Sample_ID
        refY <- coverage[sample == refSample, y]
        coverage <- coverage[sample != refSample, .(
            x = x,
            y = y - refY
        ), by = .(sample)][, .SD, .SDcols = c("x","y","sample")]    
        samplesToPlot <- samplesToPlot[2:length(samplesToPlot)]      
    }

    if(Combine_Samples){
        samplesToPlot <- list(combined = 1)
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
    image <- mdiTrackImage(layout, height, message = getBrowserTrackSetting(track, "Track_Options", "Track_Name", "svPore_coverage"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))

        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = Plot_As, #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        abline(
            h = if(Plot_As == "Read Depth") seq(0, ylim[2], 10) else -Max_Copy_Number:Max_Copy_Number, 
            col = CONSTANTS$plotlyColors$grey
        )

        I <- if(length(samplesToPlot) > 1) coverage[, sample(.N)] else TRUE
        colI <- coverage[I, unlist(samplesToPlot[sample])]
        points(coverage[I], pch = 19, col = unlist(CONSTANTS$plotlyColors[colI]))

        colI <- unlist(samplesToPlot)
        trackLegend(
            track, coord, ylim, 
            legend = names(samplesToPlot), pch = 19, col =  unlist(CONSTANTS$plotlyColors[colI])
        )
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svPore_coverageTrack <- function(track, click){
}
hover.svPore_coverageTrack <- function(track, hover){
}
brush.svPore_coverageTrack <- function(track, brush){
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_coverageTrack <- function(...) svPore_trackItems(...)
