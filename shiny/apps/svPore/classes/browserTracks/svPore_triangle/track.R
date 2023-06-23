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

    items <- track$settings$items()
    uploadName <- appStepNamesByType$upload
    samples <- as.data.table(app[[uploadName]]$outcomes$samples())

    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    ylim <- c(0, coord$width)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "Junction Clusters", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"
        for(i in seq_along(items)){
            x <- getJunctionClusters(item = items[[i]], coord = coord, rangeType = "center")

            dmsg()
            dstr(x)

            points(
                x$refPosCenter, 
                x$size,
                pch = 19,
                cex = 0.25
            )

# Classes ‘data.table’ and 'data.frame':  86602 obs. of  22 variables:
#  $ clusterN             : num  1 2 3 4 5 6 7 8 9 10 ...
#  $ edgeType             : chr  "V" "V" "D" "D" ...
#  $ eventSize            : int  117668262 8692 98319533 74228034 6259 14735 0 0 0
#  0 ...
#  $ cChromIndex1         : int  10 10 10 10 10 10 10 10 10 10 ...
#  $ cChromIndex2         : int  10 10 10 10 10 10 11 12 12 12 ...
#  $ cStrand1             : int  -1 1 1 1 1 1 1 1 1 1 ...
#  $ cStrand2             : int  1 -1 1 1 1 1 -1 -1 -1 -1 ...
#  $ cRefPos1             : int  127004201 37137222 17203466 23997033 39328322 394
# 46704 1945052 116033555 78740899 97254217 ...
#  $ cRefPos2             : int  9335939 37145914 115522999 98225067 39334581 3943
# 1969 12149997 31364182 12102774 50310792 ...
#  $ insertSize           : int  -520 -2726 3 0 -4436 -4332 0 -6 -506 0 ...       
#  $ mapQ                 : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ gapCompressedIdentity: num  0.988 0.991 0.992 0.992 0.965 ...
#  $ baseQual             : num  41.7 36.6 NA NA 32.9 37.2 NA 50 40.9 NA ...      
#  $ alnBaseQual          : num  37 35.2 35.9 38.8 33.3 ...
#  $ alnSize              : int  30236 7979 3206 24925 17617 13893 2683 6026 3002 
# 1511 ...
#  $ samples              : chr  "GM24385_HG002" "GM24385_HG002" "GM24385_HG002" "
# GM24385_HG002" ...
#  $ nSamples             : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ nInstances           : int  2 2 2 2 2 4 2 2 2 2 ...
#  $ nCanonical           : int  1 2 1 1 0 1 1 1 1 1 ...
#  $ nNonCanonical        : int  1 0 1 1 2 3 1 1 1 1 ...
#  $ GM24385_HG002        : int  2 2 2 2 2 4 2 2 2 2 ...
#  $ segmentClusterN      : num  NA 106 NA NA 108 25 NA NA NA NA ...
#  - attr(*, ".internal.selfref")=<externalptr>
#  - attr(*, "sorted")= chr "clusterN"
        }
        # 
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
