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
    jxns[, y := svx_jxnTypes[edgeType, lineN] + 0.5 - (1:.N)/.N, by = .(edgeType)]
    svx_node_plotJxnNodes(jxns[pos1In == TRUE], "pos1")
    svx_node_plotJxnNodes(jxns[pos2In == TRUE], "pos2")
    svx_node_plotJxnLines(jxns[(pos1In == TRUE | pos2In == TRUE) & edgeType != svx_edgeTypes$TRANSLOCATION], lwd, idCol)
    jxns
}
svx_node_plotGenomeJxns <- function(jxns, lwd, idCol){
    jxns[, y := size]
    svx_node_plotJxnNodes(jxns, "pos1")
    svx_node_plotJxnNodes(jxns, "pos2")
    svx_node_plotJxnLines(jxns, lwd, idCol)
    jxns
}
build.svx_nodes_track <- function(track, reference, coord, layout, trackBuffer, loadFn, idCol){
    req(coord, coord$chromosome)
    isWholeGenome <- coord$chromosome == "all"    

    # get the data to plot
    selectedSources <- getSourcesFromTrackSamples(track$settings$items())

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
            ylim = ylim, ylab = "Unique Junctions", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        jxns <- svx_getTrackJunctions(track, selectedSources, loadFn, coord, "endpoint", chromOnly = FALSE)[order(if(isWholeGenome) sample(.N) else -size)]
        if(Color_By == "sample") jxns <- dt_colorBySelectedSample(jxns, selectedSources)
        jxns[, lty := ifelse(
            edgeType == svx_edgeTypes$INVERSION, 
            ifelse(cStrand1 == 1, 1, 2), 
            1
        )]
        jxns <- if(isWholeGenome) svx_node_plotGenomeJxns(jxns, Line_Width, idCol) 
                             else svx_node_plotChromosomeJxns(jxns, Line_Width, idCol)
        svx_junctionsLegend(track, coord, ylim, selectedSources)
        trackBuffer[[track$id]] <<- jxns
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}
