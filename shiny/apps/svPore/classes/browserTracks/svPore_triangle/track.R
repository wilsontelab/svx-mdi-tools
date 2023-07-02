#----------------------------------------------------------------------
# svPore_triangle trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svPore_triangleTrackBuffer <- list()
svPore_triangleExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svPore_triangleTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click  = TRUE,
        hover  = FALSE,
        brush  = FALSE,
        items  = TRUE,
        expand = svPore_triangleExpand,
        NULL
    )
}

# build method for the S3 class; REQUIRED
build.svPore_triangleTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)
    samplesToPlot <- track$settings$items()
    sourcesToPlot <- getSvPoreSampleSources(samplesToPlot)
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    Max_SV_Size <- getBrowserTrackSetting(track, "SV_Filters", "Max_SV_Size", 0)
    if(Max_SV_Size == 0) Max_SV_Size <- as.numeric(coord$width)
    ylim <- Max_SV_Size * c(-0.05, 1.05)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = getBrowserTrackSetting(track, "Track_Options", "Track_Name", "svPore_triangle"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "Junction Clusters", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        sourceI <- 1
        jc <- do.call(rbind, lapply(names(sourcesToPlot), function(sourceId){
            jc <- applySettingsToJCs(sourceId, sourcesToPlot[[sourceId]]$Sample_ID, track) %>%
                  filterJCsByRange(coord, "center")
            jc <- cbind(
                jc[, .SD, .SDcols = c("center","size","color","cex",
                                       "nInstances","clusterN",
                                       "edgeType","cChrom1","cRefPos1","cRefPos2","node1","node2")], 
                sourceId = if(nrow(jc) == 0) character() else sourceId
            )
            sourceI <<- sourceI + 1
            jc 
        }))[sample(.N)]
        points(
            jc$center, 
            jc$size,
            pch = 19,
            cex = jc$cex,
            col = jc$color
        )        
        svPore_triangleTrackBuffer[[track$id]] <<- jc
    })

    # return the track's magick image and associated metadata
    # svPore_triangleExpand(NULL)
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svPore_triangleTrack <- function(track, click){
    jc <- svPore_triangleTrackBuffer[[track$id]]
    req(nrow(jc) > 0)    
    dist <- jc[, sqrt((center - click$coord$x) ** 2 + (size - click$coord$y) ** 2)]
    jc <- jc[which.min(dist)]    
    if(click$keys$ctrl){
        svPore_triangleExpand(jc)
        app$browser$expandingTrackId(track$id)         
    } else {
        jumpToJunctionCluster(jc)       
    }
}
hover.svPore_triangleTrack <- function(track, hover){
    # custom actions
}
brush.svPore_triangleTrack <- function(track, brush){
    # custom actions
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_triangleTrack <- function(...) svPore_trackItems(...)

# expand method for the S3 class
expand.svPore_triangleTrack <- function(track, reference, coord, layout){
    svPore_expandJunctionCluster(svPore_triangleExpand(), track, layout)
}
