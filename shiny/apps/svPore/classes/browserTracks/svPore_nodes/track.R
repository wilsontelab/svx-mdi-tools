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
svPore_nodesY <- list(

)
build.svPore_nodesTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)
    samplesToPlot <- track$settings$items()
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- c(0.5,4.5)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = getBrowserTrackSetting(track, "Track_Options", "Track_Name", "svPore_nodes"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "Junction Nodes", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        jc <- do.call(rbind, lapply(seq_along(samplesToPlot), function(sampleI){
            sample <- samplesToPlot[[sampleI]]
            sourceId <- getSourceIdFromSample_svPore(sample)
            jc <- applySettingsToJCs(sample, track) %>%
                  filterJCsByRange(coord, "endpoint", chromOnly = FALSE)
            jc[, instanceI := .I]
            jc <- cbind(
                jc[, .SD, .SDcols = c("pos1","pos2","pos1In","pos2In","edgeType","size","color","cex","instanceI","clusterN")], 
                sampleI  = if(nrow(jc) == 0) integer()   else sampleI,
                sourceId = if(nrow(jc) == 0) character() else sourceId
            )
            svPore_nodesTrackBuffer[[track$id]] <<- if(sampleI == 1) jc else rbind(svPore_nodesTrackBuffer[[track$id]], jc)
            jc 
        }))[sample(.N)]

        jc[, y := jitter(unlist(junctionTypeLines[edgeType]), amount = 0.4)]
        jc[, isIn := pos1In]
        points(
            jc[isIn == TRUE, pos1], 
            jc[isIn == TRUE, y],
            pch = 19,
            cex = jc[isIn == TRUE, cex],
            col = jc[isIn == TRUE, color]
        )

        jc[, isIn := pos2In]
        points(
            jc[isIn == TRUE, pos2], 
            jc[isIn == TRUE, y],
            pch = 19,
            cex = jc[isIn == TRUE, cex],
            col = jc[isIn == TRUE, color]
        )

        jc[, isIn := pos1In & pos2In]
        jc[isIn == TRUE, {
            lines(
                c(pos1, pos2), 
                rep(y, 2),
                col = color
            )
        }, by = .(sampleI, instanceI)]
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