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
    chrom <- coord$Chromosome
    # x <- getSampleCache(caller, cacheType, sourceId, chrom)

    x <- getSampleCache(caller, cacheType, sourceId)
    cell <- x$cells[[cellId]]
    shapeModel <- track$settings$get("Plot_Options", "Shape_Model")
    shapeKey <- tolower(shapeModel)
    replicationModel <- track$settings$get("Plot_Options", "Replication_Model")
    forceSequential <- cell$cellIsReplicating && replicationModel == "Sequential"
    repKey <- if(!cell$cellIsReplicating || forceSequential) "sequential" else "composite"
    w <- x$windows[[cell$windowPower + 1]]
    cw <- cell$windows[[shapeKey]]
    if(is.null(cw)) {
        shapeKey <- "unshaped" # this sample was not analyzed with the requested shape correction, fall back
        cw <- cell$windows[[shapeKey]]
    } 
    chromI <- w[, chrom == coord$chromosome]

    # use generic methods and any other custom code to determine the track's (dynamic) Y span
    padding <- padding(track, layout)
    height <- height(track, 3) + padding$total # or set a known, fixed height in inches
    maxCN <- 6
    ylim <- c(0, cw$RPA * maxCN)
    # ylim <- c(0, maxCN)

    pointOpacity <- 1
    defaultPointColor <- rgb(0, 0, 0, pointOpacity)
    getCnColor <- function(CN){
        cols <- c(                            # by CN (regardless of ploidy)
            defaultPointColor,                # 0 = black/grey (absence of color/copies...)
            rgb(0, 0, 1, pointOpacity),       # 1 = blue ("cool" colors are losses)
            rgb(0.1, 0.8, 0.1, pointOpacity), # 2 = subdued green ("good/go" for typical CN neutral)
            rgb(1, 0, 0, pointOpacity),       # 3 = red ("hot" colors are gains)
            rgb(1, 0.65, 0, pointOpacity),    # 4 = orange
            defaultPointColor                 # 5 = back to black/grey to make it obvious
        ) 
        col <- cols[CN + 1]
        col[is.na(col)] <- defaultPointColor
        col
    }
    minGciColor <- 30
    maxGciColor <- 55
    gcPalette <- colorRampPalette(c("orange", "blue"))(maxGciColor - minGciColor + 1)
    getGcColor <- function(gc_fraction){
        gcis <- round(gc_fraction * 100, 0)
        sapply(gcis, function(gci) {
            if(is.na(gci)) NA
            else if(gci < minGciColor) minGciColor
            else if(gci > maxGciColor) maxGciColor
            else gcPalette[gci - minGciColor + 1] 
        }) 
    }

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = "# Reads",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        points(w[chromI, start], cw$NR_wms[chromI], pch = 16, cex = 1, col = getCnColor(cw[[repKey]]$NAR[chromI] + 1))
        # points(w[chromI, start], cw[[repKey]]$CN[chromI], pch = 16, cex = 1, col = getCnColor(cw[[repKey]]$HMM[chromI] + 1))
        # points(w[chromI, start], rep(cw$RPA * 0.5, sum(chromI)), pch = 16, cex = 1, col = getGcColor(w[chromI, gc_fraction]))
        points(w[chromI, start], rep(cw$RPA * 0.7, sum(chromI)), pch = 16, cex = 1, col = getCnColor(cw[[repKey]]$HMM[chromI]))
        # points(w[chromI, start], rep(cw$RPA * 0.8, sum(chromI)), pch = 16, cex = 1, col = getCnColor(cw$NAR[chromI] + 1))

        # lines(w$start, x$hmm[[cellI]], lwd = 2, col = "red")
        abline(h = cw$RPA * 0:maxCN, col = "grey")
    })

    # return the track's magick image and associated metadata
    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )

    # cellI <- which(x$cell_id == cellId)
    # wx <- paste("w", x$window_size[cellI], sep = "_")
    # w <- x$windows[[wx]]

    # # use generic methods and any other custom code to determine the track's (dynamic) Y span
    # padding <- padding(track, layout)
    # height <- height(track, 3) + padding$total # or set a known, fixed height in inches
    # ylim <- c(0, 4)

    # # use the mdiTrackImage helper function to create the track image
    # mai <- NULL
    # image <- mdiTrackImage(layout, height, function(...){
    #     mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))

    #     # cn <- x$cn[[cellI]]
    #     # mcn <- x$windowMedians$cn[[wx]]

    #     # mcnq <- quantile(mcn, c(0.05, 0.95), na.rm = TRUE)
    #     # mcnw <- which(mcn >= mcnq[1] & mcn <= mcnq[2])
    #     # mcnf <- mcn[mcnw]
    #     # mcnf2 <- mcnf ** 2
    #     # mcnf3 <- mcnf ** 3
    #     # cnf <- cn[mcnw]
    #     # formula <- cnf ~ mcnf # + mcnf2 + mcnf3
    #     # fit <- loess(formula)

    #     # plot(0, 0, type = "n", bty = "n",
    #     #     xlim = c(0.5, 3.5), xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
    #     #     ylim = c(0.5, 3.5),  ylab = "Copy Number",
    #     #     xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

    #     # points(mcn, cn, pch = 16, cex = 0.25)
    #     # abline(0, 1, col = "blue")
    #     # points(mcnf, predict(fit, newdata = data.frame(mcnf = mcnf)), col = "red")
    #     # abline(v=mcnq)

    #     plot(0, 0, type = "n", bty = "n",
    #         xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
    #         ylim = ylim,  ylab = "Copy Number",
    #         xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"
    #     points(w$start, x$cn[[cellI]], pch = 16, cex = 1)
    #     lines(w$start, x$hmm[[cellI]], lwd = 2, col = "red")
    #     abline(h=ylim[1]:ylim[2])

    #     # plotting actions go here, e.g. points(x, y)
    # })

    # # return the track's magick image and associated metadata
    # list(
    #     ylim  = ylim,
    #     mai   = mai,
    #     image = image
    # )

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
