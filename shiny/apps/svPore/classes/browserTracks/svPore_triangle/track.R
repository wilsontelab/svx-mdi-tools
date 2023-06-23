#----------------------------------------------------------------------
# svPore_triangle trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_svPore_triangleTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE
    )
}

# build method for the S3 class; REQUIRED
build.svPore_triangleTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)
    isWholeGenome <- coord$chromosome == "all"
    samplesToPlot <- track$settings$items()

    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- as.numeric(coord$width) * c(-0.05, 1)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = track$settings$get("Track_Options","Track_Name"), function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "Junction Clusters", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # TODO: improve interim caching of a track to speed plotting on window changes    
        jc <- do.call(rbind, lapply(seq_along(samplesToPlot), function(i){
            sample <- samplesToPlot[[i]]
            jc <- loadJunctionClusters(sample = sample) %>% 
                 filterJCsBySample(sample) %>% 
                 filterJCsByRange(coord, "center") %>%
                 filterJCsBySettings(track$settings) %>%
                 svPore$setJunctionPointColors(track$settings) %>%
                 svPore$setJunctionPointSizes(track$settings)

            cbind(jc[, .SD, .SDcols = c("nodeCenter","refPosCenter","size","color","cex")], sample = sample$Sample_ID)      
        }))[sample(.N)]
        points(
            as.numeric(if(isWholeGenome) jc$nodeCenter else jc$refPosCenter), 
            jc$size,
            pch = 19,
            cex = jc$cex,
            col = jc$color
        )
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
click.svPore_triangleTrack <- function(track, x, y){
    # custom actions
}
hover.svPore_triangleTrack <- function(track, x, y){
    # custom actions
}
brush.svPore_triangleTrack <- function(track, x1, y1, x2, y2){
    # custom actions
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_triangleTrack <- function(track, session, input, reference){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select Samples",
        itemTypePlural = "Samples",
        tableData = reactive({
            uploadName <- appStepNamesByType$upload
            x <- as.data.table(app[[uploadName]]$outcomes$samples())
            req(x)
            x[, .(Sample_ID, Project)]
        }),
        keyColumn = "Sample_ID",
        extraColumns = c("Project"),
        # ,
        # options = list(
        #     XXX = list(
        #         type = "selectInput", # or textInput, etc.
        #         args = list(
        #             choices = c("aaa", "bbb"),
        #             selected = "aaa",
        #             width = "50px"                  
        #         )
        #     )
        # ),
        size = "l" # xl
    )
}
