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
                                 isMultiSample = TRUE, sampleNameFn = NULL, jxnFilterFn = NULL){
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
    Max_Y_BP    <- getBrowserTrackSetting(track, "Arcs", "Max_Y_BP", "window_width")

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svx_arcs", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        jxns <- svx_getTrackJunctions(
            track, selectedTargets, loadFn, coord, "endpoint", 
            chromOnly = FALSE, isMultiSample = isMultiSample, family = "Arcs"
        )
        if(!is.null(jxnFilterFn)) jxns <- jxnFilterFn(jxns, track) # apply app-specific filters
        if(Color_By == "sample") jxns <- dt_colorBySelectedSample(jxns, selectedTargets, isMultiSample)
        ymax <- switch(
            Max_Y_BP,
            window_width = diff(coord$range) * 0.51,
            max_sv_width = jxns[, max(abs(cRefPos1 - cRefPos2))] * 0.51,
            fixed = as.integer(getBrowserTrackSetting(track, "Arcs", "Fixed_Y_BP", 500))
        )
        ylim <<- c(0, ymax)
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = ylab(track, ""), yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 

        x <- seq(0, pi, length.out = 25)

        if(nrow(jxns) == 0) trackNoData(coord, ylim, "no matching junctions in window") else {
            jxns <- jxns[sample.int(.N, replace = FALSE)]
            jxns[, randomI := 1:.N]
            jxns[, {
                pos <- range(cRefPos1, cRefPos2)
                halfsize <- (pos[2] - pos[1] + 1) / 2
                center <- pos[1] + halfsize
                y <- sin(x) * halfsize
                if(Upside_Down) y <- ymax - y
                lines(
                    x = cos(x) * halfsize + center, 
                    y = y,
                    col = color,
                    lwd = Line_Width
                )      
            }, by = .(randomI)] # idCol
        }
        svx_junctionsLegend(track, coord, ylim, selectedTargets, isMultiSample, jxns, "Arcs", sampleNameFn)
        trackBuffer[[track$id]] <- jxns
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}
