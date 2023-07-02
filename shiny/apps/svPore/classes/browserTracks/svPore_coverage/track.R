#----------------------------------------------------------------------
# svPore_nodes trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svPore_nodesExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svPore_nodesTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = TRUE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        expand = svPore_nodesExpand,
        NULL
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
    jc
}
plotGenomeNodes <- function(jc){
    jc[, y := size]
    plotSvNodeEndpoints(jc, "pos1")
    plotSvNodeEndpoints(jc, "pos2")
    plotSvNodeLines(jc)
    jc
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
                jc[, .SD, .SDcols = c("pos1","pos2","pos1In","pos2In","size","color","cex","clusterN",
                                      "edgeType","cChrom1","cRefPos1","cRefPos2","node1","node2")], 
                sourceId = if(nrow(jc) == 0) character() else sourceId
            )
            sourceI <<- sourceI + 1
            jc 
        }))[order(if(isWholeGenome) sample(.N) else -size)]
        jc <- if(isWholeGenome) plotGenomeNodes(jc) else plotChromosomeNodes(jc)
        svPore_nodesTrackBuffer[[track$id]] <<- jc
    })

    # return the track's magick image and associated metadata
    # svPore_nodesExpand(NULL)
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svPore_nodesTrack <- function(track, click){
    jc <- svPore_nodesTrackBuffer[[track$id]]
    req(nrow(jc) > 0)  
    y2 <- jc[, ((y - click$coord$y) / click$coord$y) ** 2]
    dist <- jc[, pmin(
        sqrt(((pos1 - click$coord$x) / click$coord$x) ** 2 + y2),
        sqrt(((pos2 - click$coord$x) / click$coord$x) ** 2 + y2)
    )]
    jc <- jc[which.min(dist)]
    if(click$keys$ctrl){
        svPore_nodesExpand(jc)
        app$browser$expandingTrackId(track$id)         
    } else {
        jumpToJunctionCluster(jc)       
    }
}
hover.svPore_nodesTrack <- function(track, hover){
    # custom actions
}
brush.svPore_nodesTrack <- function(track, brush){
    # custom actions
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_nodesTrack <- function(...) svPore_trackItems(...)

# expand methods for the S3 class
expand.svPore_nodesTrack <- function(track, reference, coord, layout){
    svPore_expandJunctionCluster(svPore_nodesExpand(), track, layout)
}
