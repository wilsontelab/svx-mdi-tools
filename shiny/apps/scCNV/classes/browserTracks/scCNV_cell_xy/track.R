#----------------------------------------------------------------------
# scCNV_cell_xy trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
# CONSTANTS$scCNV$binSize <- 2e4
# input.scCNV_cell_xyTrack <- function(track, name) paste("scCNV_cell_xy", name, track$id, sep = "_")

# constructor for the S3 class
new_scCNV_cell_xyTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        genome = FALSE, # whether the track supports a whole genome view
        adjustsWidth = FALSE
    )
}

# build method for the S3 class; REQUIRED
build.scCNV_cell_xyTrack <- function(track, reference, coord, layout){
    cellId <- track$settings$get('Plot_Options', 'Cell_ID') 
    req(cellId)

    caller <- paste("scCNV_heatmapTrack", track$id, sep = "-")
    cacheType <- "working"
    sourceId <- track$settings$get('Plot_Options', 'Sample')    
    chrom <- coord$chromosome    
    x <- getSampleCache(caller, cacheType, sourceId, chrom)


    cellI <- which(x$cell_id == cellId)
    w <- x$windows[[paste("w", x$window_size[cellI], sep = "_")]]


    # use generic methods and any other custom code to determine the track's (dynamic) Y span
    padding <- padding(track, layout)
    height <- height(track, 1) + padding$total # or set a known, fixed height in inches
    ylim <- c(-2, 2)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "Copy Number",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        dprint(coord$range)
        dprint(range(w$start))
        dprint(ylim)
        dprint(range(x$cn[[cellI]]))
        points(w$start, x$cn[[cellI]], pch = 16, cex = 1)
        lines(w$start, x$hmm[[cellI]], lwd = 2, col = "red")
        abline(h=ylim[1]:ylim[2])
            

        # plotting actions go here, e.g. points(x, y)
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )

# List of 5
#  $ cell_id    : chr [1:10] "0" "1" "2" "3" ...
#  $ modal_CN   : int [1:10] 2 2 2 2 2 2 2 2 2 2
#  $ rejected   : logi [1:10] TRUE TRUE FALSE FALSE FALSE TRUE ...
#  $ window_size: int [1:10] 1 1 3 3 3 1 1 3 3 3
#  $ z          :List of 2
#   ..$ cn :List of 3
#   .. ..$ : num [1:12446, 1:10] NA NA NA NA NA NA NA NA NA NA ...
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : NULL
#   .. .. .. ..$ : chr [1:10] "0" "1" "2" "3" ...
#   .. ..$ : num [1:12446, 1:10] NA NA NA NA NA NA NA NA NA NA ...
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : NULL
#   .. .. .. ..$ : chr [1:10] "0" "1" "2" "3" ...
#   .. ..$ : num [1:12446, 1:10] NA NA NA NA NA NA NA NA NA NA ...
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : NULL
#   .. .. .. ..$ : chr [1:10] "0" "1" "2" "3" ...
#   ..$ hmm:List of 3
#   .. ..$ : num [1:12446, 1:10] NA NA NA NA NA NA NA NA NA NA ...
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : NULL
#   .. .. .. ..$ : chr [1:10] "0" "1" "2" "3" ...
#   .. ..$ : num [1:12446, 1:10] NA NA NA NA NA NA NA NA NA NA ...
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : NULL
#   .. .. .. ..$ : chr [1:10] "0" "1" "2" "3" ...
#   .. ..$ : num [1:12446, 1:10] NA NA NA NA NA NA NA NA NA NA ...
#   .. .. ..- attr(*, "dimnames")=List of 2
#   .. .. .. ..$ : NULL
#   .. .. .. ..$ : chr [1:10] "0" "1" "2" "3" ...

    # req(d)   
    # caller <- paste("scCNV_cell_xyTrack", track$id, sep = "-")
    # cacheType <- "working"
    # sourceId <- track$settings$get('Plot_Options', 'Sample')    
    # chrom <- coord$chromosome    
    # x <- getSampleCache(caller, cacheType, sourceId, chrom)

    # keptCells <- which(!x$rejected) # TODO: filter cells should obey user overrides
    # nCells <- length(keptCells)
    # maxZ <- 2 # TODO: expose as setting?

    # bins <- d$minBinN:min(d$maxBinN, d$maxBin)
    # isOverflow <- d$maxBinN > d$maxBin
    # if(isOverflow) padding <- matrix(NA, nrow = d$maxBinN - d$maxBin, ncol = length(keptCells))
    # plotType <- track$browser[[ input.scCNV_cell_xyTrack(track, "plotType") ]]
    # getVal <- function(cnc){
    #     x <- x$z[[plotType]][[cnc + 2]][bins, keptCells]
    #     if(isOverflow) x <- rbind(x, padding)
    #     x
    # }

    # cellSort <- track$browser[[ input.scCNV_cell_xyTrack(track, "cellSort") ]]
    # keptCells <- switch(
    #     cellSort,
    #     resolution = keptCells[order(-x$window_size[keptCells], keptCells)],
    #     hclust = { # TODO: expose methods as settings  
    #         z <- x$z$cn[[0 + 2]][bins, keptCells]
    #         # z <- z[apply(z, 1, function(v) !all(is.na(v))), ]
    #         dst <- dist(t(z), method = "euclidean")
    #         keptCells[hclust(dst, method = "complete")$order]
    #     },
    #     keptCells
    # )

    # fileName <- paste0(app$NAME, "-", track$id, "-scCNV_cell_xyTrack.png")     
    # pngFile <- file.path(sessionDirectory, fileName)
    # yScaleFactor <- track$settings$get("Plot_Options", "Cell_Height_Pixels")
    # yScaleFactor <- if(is.null(yScaleFactor) || yScaleFactor == "auto"){
    #     minHeight <- track$settings$get("Plot_Options", "Min_Height_Pixels")
    #     if(nCells > minHeight || nCells == 0) 1 else ceiling(minHeight / nCells)
    # } else as.integer(yScaleFactor)
    # if(plotType == "cn"){
    #     loss    <- zToIntensity_neg_decay(-getVal(-1), maxZ = maxZ)
    #     neutral <- zToIntensity_symmetric( getVal( 0), maxZ = maxZ)
    #     gain    <- zToIntensity_neg_decay( getVal( 1), maxZ = maxZ)
    # } else {
    #     loss    <- getVal(-1)
    #     neutral <- getVal( 0)
    #     gain    <- getVal( 1)
    # }
    # saveHeatMap_three_color(
    #     loss,
    #     neutral,
    #     gain,
    #     file = pngFile, 
    #     mirror = FALSE, 
    #     xScaleFactor = 1 / d$binsPerPixel, 
    #     yScaleFactor = yScaleFactor # TODO: either auto or hard number from settings
    # ) 
    # list(image = pngToMdiTrackImage(pngFile, layout))
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.scCNV_cell_xyTrack <- function(track, x, y){
    # custom actions
}
hover.scCNV_cell_xyTrack <- function(track, x, y){
    # custom actions
}
brush.scCNV_cell_xyTrack <- function(track, x1, y1, x2, y2){
    # custom actions
}

# # method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# # used when a track can take a list of items to be plotted together and the item list icon is clicked
# items.scCNV_cell_xyTrack <- function(track, session, input, reference){
#     showTrackItemsDialog(
#         settings,
#         session,
#         title = "Select XYZ",
#         itemTypePlural = "Items",
#         tableData = function() data.frame(xxx = 1, yyy = 2),
#         keyColumn = "xxx",
#         extraColumns = c("yyy"),
#         options = list(
#             XXX = list(
#                 type = "selectInput", # or textInput, etc.
#                 args = list(
#                     choices = c("aaa", "bbb"),
#                     selected = "aaa",
#                     width = "50px"                  
#                 )
#             )
#         ),
#         size = "l"
#     )
# }

# # implement track-level navigation
# observers.scCNV_cell_xyTrack <- list()
# navigation.scCNV_cell_xyTrack <- function(track, session){

#     # # handle nav observers, e.g., buttons
#     # for(x in observers.scCNV_cell_xyTrack) x$destroy()
#     # input <- session$input
#     # observers.scCNV_cell_xyTrack$prevPage <<- observeEvent(input$prevPage, {
#     #     pageNumber <- as.integer(input$pageNumber)
#     #     if(pageNumber == 1) return()
#     #     updateTextInput(session, "pageNumber", value = pageNumber - 1)
#     # })
#     # observers.scCNV_cell_xyTrack$nextPage <<- observeEvent(input$nextPage, {
#     #     pageNumber <- as.integer(input$pageNumber)
#     #     # if(pageNumber == 1) return()
#     #     updateTextInput(session, "pageNumber", value = pageNumber + 1)
#     # })

#     # return the navigation UI elements
#     tagList(
#         tags$div(
#             class = "trackBrowserInput",
#             radioButtons(
#                 session$ns(input.scCNV_cell_xyTrack(track, "plotType")),
#                 "Plot Type",
#                 choices = c("cn", "hmm"),
#                 selected = "cn",
#                 inline = TRUE
#             )
#         ),
#         tags$div(
#             class = "trackBrowserInput",
#             style = "margin-left: 20px;",
#             radioButtons(
#                 session$ns(input.scCNV_cell_xyTrack(track, "cellSort")),
#                 "Cell Sort",
#                 choices = c("resolution", "hclust"),
#                 selected = "resolution",
#                 inline = TRUE
#             )
#         )
#     )
# }
