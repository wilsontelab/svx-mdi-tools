#----------------------------------------------------------------------
# scCNV_heatmap trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
data.scCNV_heatmapTrack <- list()
CONSTANTS$scCNV$binSize <- 2e4
input.scCNV_heatmapTrack <- function(track, name) paste("scCNV_heatmap", name, track$id, sep = "_")

# constructor for the S3 class
new_scCNV_heatmapTrack <- function(trackId) {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = TRUE, # whether the track offers a custom, additional row of within-track navigation inputs
        genome = FALSE, # whether the track supports a whole genome view
        adjustsWidth = TRUE
    )
}
delete.scCNV_heatmapTrack <- function(track){
    # caller <- paste("scCNV_heatmapTrack", track$id, sep = "-")
    # purgeSampleCache(samples.scCNV_heatmapTrack[[caller]], caller)
}

# pre-plot action to override browser coordinates and width
adjustWidth.scCNV_heatmapTrack <- function(track, reference, coord, layout){
    data.scCNV_heatmapTrack[[track$id]] <<- NULL
    nullResult <- list(coord = coord, layout = layout)
    caller <- paste("scCNV_heatmapTrack", track$id, sep = "-")
    cacheType <- "common"
    sourceId <- track$settings$get('Plot_Options', 'Sample')
    if(is.null(sourceId) || is.na(sourceId) || sourceId == "") return(nullResult)
    x <- getSampleCache(caller, cacheType, sourceId)
    isAll <- coord$chromosome == "all"
    binN <- if(isAll) x$rowRanges[cumsum(as.numeric(start)) <= coord$end & cumsum(as.numeric(end)) >= coord$start, bin_n]
                 else x$rowRanges[chrom == coord$chromosome & start <= coord$end & end >= coord$start, chrom_bin_n]

    if(length(binN) == 0) return(nullResult)

    nBins <- diff(range(binN)) + 1
    maxWidth <- layout$plotWidth * layout$dpi
    binsPerPixel <- ceiling(nBins / maxWidth)
    width <- ceiling(nBins / binsPerPixel)
    nBins <- binsPerPixel * width
    if(binsPerPixel == 1) while(width + nBins <= maxWidth) width <- width + nBins
    minBinN <- binN[1]
    maxBinN <- minBinN + nBins - 1
    startPos <- (minBinN - 1) * CONSTANTS$scCNV$binSize + 1
    endPos <- startPos + CONSTANTS$scCNV$binSize * nBins - 1
    data.scCNV_heatmapTrack[[track$id]] <<- list(
        maxBin = max(binN),
        minBinN = minBinN,
        maxBinN = maxBinN,
        nBins = nBins,
        width = width,
        binsPerPixel = nBins / width,
        startPos = startPos,
        endPos = endPos
    )
    
    list(
        coord = getBrowserCoord(coord$chromosome, startPos, endPos),
        layout = adjustLayoutWidth(layout, width)
    )
}

# build method for the S3 class; REQUIRED
build.scCNV_heatmapTrack <- function(track, reference, coord, layout){
    d <- data.scCNV_heatmapTrack[[track$id]]
    req(d)   
    caller <- paste("scCNV_heatmapTrack", track$id, sep = "-")
    cacheType <- "working"
    sourceId <- track$settings$get('Plot_Options', 'Sample')    
    chrom <- coord$chromosome    
    x <- getSampleCache(caller, cacheType, sourceId, chrom)

    keptCells <- which(!x$rejected) # TODO: filter cells should obey user overrides
    nCells <- length(keptCells)
    maxZ <- 2 # TODO: expose as setting?

    bins <- d$minBinN:min(d$maxBinN, d$maxBin)
    isOverflow <- d$maxBinN > d$maxBin
    if(isOverflow) padding <- matrix(NA, nrow = d$maxBinN - d$maxBin, ncol = length(keptCells))
    plotType <- track$browser[[ input.scCNV_heatmapTrack(track, "plotType") ]]
    getVal <- function(cnc){
        x <- x$z[[plotType]][[cnc + 2]][bins, keptCells]
        if(isOverflow) x <- rbind(x, padding)
        x
    }

    cellSort <- track$browser[[ input.scCNV_heatmapTrack(track, "cellSort") ]]
    keptCells <- switch(
        cellSort,
        resolution = keptCells[order(-x$window_size[keptCells], keptCells)],
        hclust = { # TODO: expose methods as settings  
            z <- x$z$cn[[0 + 2]][bins, keptCells]
            # z <- z[apply(z, 1, function(v) !all(is.na(v))), ]
            dst <- dist(t(z), method = "euclidean")
            keptCells[hclust(dst, method = "complete")$order]
        },
        keptCells
    )

    fileName <- paste0(app$NAME, "-", track$id, "-scCNV_heatmapTrack.png")     
    pngFile <- file.path(sessionDirectory, fileName)
    yScaleFactor <- track$settings$get("Plot_Options", "Cell_Height_Pixels")
    yScaleFactor <- if(is.null(yScaleFactor) || yScaleFactor == "auto"){
        minHeight <- track$settings$get("Plot_Options", "Min_Height_Pixels")
        if(nCells > minHeight || nCells == 0) 1 else ceiling(minHeight / nCells)
    } else as.integer(yScaleFactor)
    if(plotType == "cn"){
        loss    <- zToIntensity_neg_decay(-getVal(-1), maxZ = maxZ)
        neutral <- zToIntensity_symmetric( getVal( 0), maxZ = maxZ)
        gain    <- zToIntensity_neg_decay( getVal( 1), maxZ = maxZ)
    } else {
        loss    <- getVal(-1)
        neutral <- getVal( 0)
        gain    <- getVal( 1)
    }
    saveHeatMap_three_color(
        loss,
        neutral,
        gain,
        file = pngFile, 
        mirror = FALSE, 
        xScaleFactor = 1 / d$binsPerPixel, 
        yScaleFactor = yScaleFactor # TODO: either auto or hard number from settings
    ) 
    list(image = pngToMdiTrackImage(pngFile, layout))
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.scCNV_heatmapTrack <- function(track, x, y){
    # custom actions
}
hover.scCNV_heatmapTrack <- function(track, x, y){
    # custom actions
}
brush.scCNV_heatmapTrack <- function(track, x1, y1, x2, y2){
    # custom actions
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.scCNV_heatmapTrack <- function(track, session, input, reference){
    showTrackItemsDialog(
        settings,
        session,
        title = "Select XYZ",
        itemTypePlural = "Items",
        tableData = function() data.frame(xxx = 1, yyy = 2),
        keyColumn = "xxx",
        extraColumns = c("yyy"),
        options = list(
            XXX = list(
                type = "selectInput", # or textInput, etc.
                args = list(
                    choices = c("aaa", "bbb"),
                    selected = "aaa",
                    width = "50px"                  
                )
            )
        ),
        size = "l"
    )
}

# implement track-level navigation
observers.scCNV_heatmapTrack <- list()
navigation.scCNV_heatmapTrack <- function(track, session){

    # # handle nav observers, e.g., buttons
    # for(x in observers.scCNV_heatmapTrack) x$destroy()
    # input <- session$input
    # observers.scCNV_heatmapTrack$prevPage <<- observeEvent(input$prevPage, {
    #     pageNumber <- as.integer(input$pageNumber)
    #     if(pageNumber == 1) return()
    #     updateTextInput(session, "pageNumber", value = pageNumber - 1)
    # })
    # observers.scCNV_heatmapTrack$nextPage <<- observeEvent(input$nextPage, {
    #     pageNumber <- as.integer(input$pageNumber)
    #     # if(pageNumber == 1) return()
    #     updateTextInput(session, "pageNumber", value = pageNumber + 1)
    # })

    # return the navigation UI elements
    tagList(
        tags$div(
            class = "trackBrowserInput",
            radioButtons(
                session$ns(input.scCNV_heatmapTrack(track, "plotType")),
                "Plot Type",
                choices = c("cn", "hmm"),
                selected = "cn",
                inline = TRUE
            )
        ),
        tags$div(
            class = "trackBrowserInput",
            style = "margin-left: 20px;",
            radioButtons(
                session$ns(input.scCNV_heatmapTrack(track, "cellSort")),
                "Cell Sort",
                choices = c("resolution", "hclust"),
                selected = "resolution",
                inline = TRUE
            )
        )
    )
}
