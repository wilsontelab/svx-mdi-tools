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
}
