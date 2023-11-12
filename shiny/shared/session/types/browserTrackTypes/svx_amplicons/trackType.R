#----------------------------------------------------------------------
# svx_amplicons is a track type for showing the spans of one or more amplicons
# using the same concepts as genome_spans tracks
#----------------------------------------------------------------------

# build method for the S3 class; REQUIRED
build.svx_amplicons_track <- function(track, reference, coord, layout, 
                                     trackBuffer, ampliconsFamily = "Amplicons",
                                     isMultiSample = FALSE){
    req(coord, coord$chromosome, coord$chromosome != "all")

    # get the data to plot
    amplicons <- track$settings$items()
    if(isMultiSample) amplicons <- getSourcesFromTrackSamples(amplicons)
    itemData <- svx_getTrackAmpliconTargets(track, amplicons = amplicons)
    itemData[, ":="(source = ampliconTargetKey, strand = ".", x1 = start, x2 = end)]
    itemData <- itemData[order(x1, x2)]
    nItems <- nrow(itemData)
    itemNames <- itemData$ampliconTargetKey

    # set the plot layout and axes
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches  
    Color_By    <- getBrowserTrackSetting(track, "Arcs", "Color_By", "edgeType")
    Line_Width  <- getBrowserTrackSetting(track, "Arcs", "Line_Width", 1)

    # convert to span lines
    itemsList <- lapply(itemNames, list)
    names(itemsList) <- itemNames
    Plot_Spans_As <- trimws(getTrackSetting(track, ampliconsFamily, "Plot_Spans_As", "packed_spans"))
    x <- if(Plot_Spans_As == "packed_spans"){
        getPackedSpans(track, coord, ampliconsFamily, itemsList, itemData)
    } else { # unpacked_spans
        getUnpackedSpans(itemsList, itemData)
    } 
    itemsList <- x$itemsList
    itemData <- x$itemData
    if(!is.null(trackBuffer)) trackBuffer[[track$id]] <- itemData
    buildSpanTrackImage (
        track, coord, layout,
        itemsList, itemNames, itemData,
        stranded = FALSE, allowNeg = FALSE, ylab = " ", ylim = c(itemsList$ymin, itemsList$ymax), yaxt = "n",
        dataFamily = ampliconsFamily, yAxisFamily = ampliconsFamily, hLines = FALSE
    )

    # # use the mdiTrackImage helper function to create the track image
    # mai <- NULL
    # image <- mdiTrackImage(layout, height, message = "svx_amplicons", function(...){
    #     mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
    #     jxns <- svx_getTrackJunctions(track, selectedTargets, loadFn, coord, "endpoint", 
    #                                   chromOnly = FALSE, isMultiSample = isMultiSample, family = "Arcs")
    #     if(!is.null(jxnFilterFn)) jxns <- jxnFilterFn(jxns, track) # apply app-specific filters
    #     if(Color_By == "sample") jxns <- dt_colorBySelectedSample(jxns, selectedTargets, isMultiSample)
    #     ymax <- switch(
    #         Max_Y_BP,
    #         window_width = diff(coord$range) * 0.51,
    #         max_sv_width = jxns[, max(abs(cRefPos1 - cRefPos2))] * 0.51,
    #         fixed = as.integer(getBrowserTrackSetting(track, "Arcs", "Fixed_Y_BP", 500))
    #     )
    #     ylim <<- c(0, ymax)
    #     plot(0, 0, type = "n", bty = "n",
    #         xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
    #         ylim = ylim, ylab = ylab(track, ""), yaxt = "n",
    #         xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 

    #     x <- seq(0, pi, length.out = 25)

    #     if(nrow(jxns) == 0) trackNoData(coord, ylim, "no matching junctions in window") else {
    #         jxns <- jxns[sample.int(.N, replace = FALSE)]
    #         jxns[, randomI := 1:.N]
    #         jxns[, {
    #             pos <- range(cRefPos1, cRefPos2)
    #             halfsize <- (pos[2] - pos[1] + 1) / 2
    #             center <- pos[1] + halfsize
    #             y <- sin(x) * halfsize
    #             if(Upside_Down) y <- ymax - y
    #             lines(
    #                 x = cos(x) * halfsize + center, 
    #                 y = y,
    #                 col = color,
    #                 lwd = Line_Width
    #             )      
    #         }, by = .(randomI)] # idCol
    #     }
    #     svx_junctionsLegend(track, coord, ylim, selectedTargets, isMultiSample, jxns, "Arcs", sampleNameFn)
    #     trackBuffer[[track$id]] <- jxns
    # })

    # # return the track's magick image and associated metadata
    # list(ylim = ylim, mai = mai, image = image)
}
