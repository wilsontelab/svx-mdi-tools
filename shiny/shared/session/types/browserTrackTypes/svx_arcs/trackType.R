#----------------------------------------------------------------------
# svx_arcs is a track type for plotting SV node pairs as semicircle "jumps" between them
# expects:
#   one or more track$settings$items from showTrackSamplesDialog()
#   dataFn(track, sourceId, samples) must return data.table(samples, pos1 = integer64, pos2 = integer64, size = integer64, cex, color, ...)
# data may arise from one or more sources, derived from Adjust Sample Sources
# each source may have one or more samples, as a result of multi-sample find actions in the relevant svX pipelines
#----------------------------------------------------------------------

# constructor
new_svx_arcsTrack <- function(trackId, expandReactive) {
    list(
        click  = TRUE,
        hover  = FALSE,
        brush  = FALSE,
        items  = TRUE,
        expand = expandReactive,
        expand2 = TRUE,
        navigation = TRUE
    )
}

# build method for the S3 class; REQUIRED
build.svx_arcs_track <- function(track, reference, coord, layout, trackBuffer, loadFn, idCol, 
                                 isMultiSample = TRUE, sampleNameFn = NULL){
    req(coord, coord$chromosome)
    isWholeGenome <- coord$chromosome == "all"    

    # get the data to plot
    selectedTargets <- track$settings$items()
    if(isMultiSample) selectedTargets <- getSourcesFromTrackSamples(selectedTargets)

    # set the plot layout and axes
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    Upside_Down <- getBrowserTrackSetting(track, "Arcs", "Upside_Down", FALSE)    
    Color_By    <- getBrowserTrackSetting(track, "Arcs", "Color_By", "edgeType")
    Line_Width  <- getBrowserTrackSetting(track, "Arcs", "Line_Width", 1)
    Max_SV_Size <- getBrowserTrackSetting(track, "Filters", "Max_SV_Size", 0)
    ymax <- (if(Max_SV_Size == 0) diff(coord$range) else Max_SV_Size) / 2
    ylim <- c(0, ymax)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svx_arcs", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "", yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        jxns <- svx_getTrackJunctions(track, selectedTargets, loadFn, coord, "endpoint", 
                                      chromOnly = FALSE, isMultiSample = isMultiSample, family = "Arcs")
        if(Color_By == "sample") jxns <- dt_colorBySelectedSample(jxns, selectedTargets, isMultiSample)
        x <- seq(0, pi, length.out = 25)
        if(nrow(jxns) == 0) trackNoData(coord, ylim, "no matching junctions in window") else {
            jxns[, {
                halfsize <- (cRefPos2 - cRefPos1 + 1) / 2
                center <- cRefPos1 + halfsize
                y <- sin(x) * halfsize
                if(Upside_Down) y <- ymax - y
                lines(
                    x = cos(x) * halfsize + center, 
                    y = y,
                    col = color,
                    # lwd = lwd
                )      
            }, by = idCol]
        }
        svx_junctionsLegend(track, coord, ylim, selectedTargets, isMultiSample, jxns, "Arcs", sampleNameFn)
        trackBuffer[[track$id]] <<- jxns
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}
