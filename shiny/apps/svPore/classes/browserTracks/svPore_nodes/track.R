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
svPore_coverageTrackBuffer <- list()

build.svPore_coverageTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)
    Max_Points   <- getBrowserTrackSetting(track, "Coverage", "Max_Points", 1000)
    Max_Coverage <- getBrowserTrackSetting(track, "Coverage", "Max_Coverage", 0) # 0 means "auto"
    if(Max_Coverage == 0) Max_Coverage <- 1e8
    samplesToPlot <- track$settings$items()
    sourcesToPlot <- getSvPoreSampleSources(samplesToPlot)
    coverage <- do.call(rbind, lapply(names(sourcesToPlot), function(sourceId){
        do.call(rbind, lapply(sourcesToPlot[[sourceId]]$Sample_ID, function(sample_){
            coverage <- filterCoverageByRange(sourceId, sample_, coord) %>% 
                        aggregateCoverageBins(Max_Points)
            coverage[, sample := sample_]
            coverage
        }))
    }))   
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- c(coverage[, min(coverage)], coverage[, min(max(coverage), Max_Coverage)]) 

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = getBrowserTrackSetting(track, "Track_Options", "Track_Name", "svPore_coverage"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "Coverage", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        for(sample_ in unique(coverage$sample)){
            x <- coverage[sample == sample_]
            points(if(coord$chromosome == "all") as.numeric(x$genomeStart) else x$start, x$coverage, pch = 19)
        }
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
