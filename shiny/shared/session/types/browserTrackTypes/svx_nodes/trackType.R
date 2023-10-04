#----------------------------------------------------------------------
# svx_nodes is a track type for plotting SV node pairs as connect lines on a horizontal map
# expects:
#   one or more track$settings$items from showTrackSamplesDialog()
#   dataFn(track, sourceId, samples) must return data.table(samples, pos1 = integer64, pos2 = integer64, size = integer64, cex, color, ...)
# data may arise from one or more sources, derived from Adjust Sample Sources
# each source may have one or more samples, as a result of multi-sample find actions in the relevant svX pipelines
#----------------------------------------------------------------------

# constructor
new_svx_nodesTrack <- function(trackId, expandReactive) {
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
svx_node_plotJxnNodes <- function(jxns, pos){
    points(
        jxns[[pos]], 
        jxns[, y],
        pch = 19,
        cex = jxns[, cex],
        col = jxns[, color]
    )       
}
svx_node_plotJxnLines <- function(jxns, lwd, idCol){
    jxns[, {
        lines(
            c(pos1, pos2), 
            rep(y, 2),
            col = color,
            lwd = lwd,
            lty = lty
        )
    }, by = idCol]
}
svx_node_plotChromosomeJxns <- function(jxns, lwd, idCol){
    jxns[, y := svx_jxnType_codeToX(edgeType, "lineN") + 0.5 - (1:.N)/.N, by = .(edgeType)]
    svx_node_plotJxnNodes(jxns[pos1In == TRUE], "pos1")
    svx_node_plotJxnNodes(jxns[pos2In == TRUE], "pos2")
    svx_node_plotJxnLines(jxns[(pos1In == TRUE | pos2In == TRUE) & edgeType != svx_edgeTypes$TRANSLOCATION], lwd, idCol)
    axisCodes <- c("D","U","V","T")
    axis(
        side = 2, 
        at = svx_jxnType_codeToX(axisCodes, "lineN"), 
        labels = svx_jxnType_codeToX(axisCodes, "name"), 
        las = 1
    )
    jxns
}
svx_node_plotGenomeJxns <- function(jxns, lwd, idCol){
    jxns[, y := size]
    svx_node_plotJxnNodes(jxns, "pos1")
    svx_node_plotJxnNodes(jxns, "pos2")
    svx_node_plotJxnLines(jxns, lwd, idCol)
    jxns
}
build.svx_nodes_track <- function(track, reference, coord, layout, trackBuffer, loadFn, idCol, 
                                  isMultiSample = TRUE, sampleNameFn = NULL, jxnFilterFn = NULL){
    req(coord, coord$chromosome)
    isWholeGenome <- coord$chromosome == "all"    

    # get the data to plot
    selectedTargets <- track$settings$items()
    if(isMultiSample) selectedTargets <- getSourcesFromTrackSamples(selectedTargets)

    # set the plot layout and axes
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    Color_By <- getBrowserTrackSetting(track, "Points", "Color_By", "edgeType")
    Line_Width <- getBrowserTrackSetting(track, "Points", "Line_Width", 1)
    ylim <- if(isWholeGenome) {
        Max_SV_Size <- getBrowserTrackSetting(track, "Filters", "Max_SV_Size", 0)
        if(Max_SV_Size == 0) coord$range else c(1, as.numeric(Max_SV_Size))
    } else c(0.45, 4.55)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svx_nodes", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = if(isWholeGenome) "SV Size" else "Unique Junctions", yaxt = if(isWholeGenome) "s" else "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        jxns <- svx_getTrackJunctions(
            track, selectedTargets, loadFn, coord, "endpoint", chromOnly = FALSE, isMultiSample = isMultiSample
        )[order(if(isWholeGenome) sample(.N) else -size)]
        if(!is.null(jxnFilterFn)) jxns <- jxnFilterFn(jxns, track) # apply app-specific filters
        if(Color_By == "sample") jxns <- dt_colorBySelectedSample(jxns, selectedTargets, isMultiSample)
        jxns[, lty := ifelse(
            edgeType == svx_edgeTypes$INVERSION, 
            ifelse(cStrand1 == 1, 1, 2), 
            1
        )]
        jxns <- if(isWholeGenome) svx_node_plotGenomeJxns(jxns, Line_Width, idCol) 
                             else svx_node_plotChromosomeJxns(jxns, Line_Width, idCol)
        if(nrow(jxns) == 0) trackNoData(coord, ylim, "no matching junctions in window")
        svx_junctionsLegend(track, coord, ylim, selectedTargets, isMultiSample, jxns, sampleNameFn = sampleNameFn)
        trackBuffer[[track$id]] <- jxns
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}
