#----------------------------------------------------------------------
# svPore_nodes trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_svPore_nodesTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE
    )
}

# build method for the S3 class; REQUIRED
svPore_nodesTrackBuffer <- list()
plotSvNodeEndpoints <- function(jc, pos){
    points(
        jc[[pos]], 
        jc[, y],
        pch = 19,
        cex = jc[, cex],
        col = jc[, color]
    )       
}
plotSvNodeLines <- function(jc){
    jc[, {
        lines(
            c(pos1, pos2), 
            rep(y, 2),
            col = color
        )
    }, by = .(clusterN)]     
}
plotChromosomeNodes <- function(jc){
    jc[, y := junctionTypeLines[[edgeType]] - 0.5 + (1:.N)/.N, by = .(edgeType)]
    plotSvNodeEndpoints(jc[pos1In == TRUE], "pos1")
    plotSvNodeEndpoints(jc[pos2In == TRUE], "pos2")
    plotSvNodeLines(jc[pos1In == TRUE & pos2In == TRUE])
}
plotGenomeNodes <- function(jc){
    jc[, y := size]
    plotSvNodeEndpoints(jc, "pos1")
    plotSvNodeEndpoints(jc, "pos2")
    plotSvNodeLines(jc)
}
build.svPore_nodesTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)
    isWholeGenome <- coord$chromosome == "all"
    sourcesToPlot <- getSvPoreSampleSources(track$settings$items())
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- if(isWholeGenome) {
        Max_SV_Size <- getBrowserTrackSetting(track, "SV_Filters", "Max_SV_Size", 0)
        if(Max_SV_Size == 0) coord$range else c(1, as.numeric(Max_SV_Size))
    } else c(0.45,4.55)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = getBrowserTrackSetting(track, "Track_Options", "Track_Name", "svPore_nodes"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "Junction Nodes", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        sourceI <- 1
        jc <- do.call(rbind, lapply(names(sourcesToPlot), function(sourceId){
            jc <- applySettingsToJCs(sourceId, sourcesToPlot[[sourceId]]$Sample_ID, track) %>%
                  filterJCsByRange(coord, "endpoint", chromOnly = FALSE)
            jc <- cbind(
                jc[, .SD, .SDcols = c("pos1","pos2","pos1In","pos2In","edgeType","size","color","cex","clusterN")], 
                sourceId = if(nrow(jc) == 0) character() else sourceId
            )
            svPore_nodesTrackBuffer[[track$id]] <<- if(sourceI == 1) jc else rbind(svPore_nodesTrackBuffer[[track$id]], jc)
            sourceI <<- sourceI + 1
            jc 
        }))[order(if(isWholeGenome) sample(.N) else -size)]
        if(isWholeGenome) plotGenomeNodes(jc) else plotChromosomeNodes(jc)
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
click.svPore_nodesTrack <- function(track, x, y){
    # jc <- svPore_nodesTrackBuffer[[track$id]]
    # dist <- jc[, sqrt((center - x) ** 2 + (size - y) ** 2)]
    # jc <- jc[which.min(dist)]
    # req(nrow(jc) > 0)
    
    # dprint(jc)
}
hover.svPore_nodesTrack <- function(track, x, y){
    # custom actions
}
brush.svPore_nodesTrack <- function(track, x1, y1, x2, y2){
    # custom actions
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_nodesTrack <- function(...) svPore_trackItems(...)
