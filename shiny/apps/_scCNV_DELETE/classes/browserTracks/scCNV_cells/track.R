#----------------------------------------------------------------------
# scCNV_cells trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_scCNV_cellsTrack <- function() {
    list( # whether the track type has `click`, `hover`, and/or `items` methods
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = FALSE,
        navigation = TRUE, # whether the track offers a custom, additional row of within-track navigation inputs
        genome = FALSE # whether the track supports a whole genome view
    )
}

# build method for the S3 class; REQUIRED
build.scCNV_cellsTrack <- function(settings, input, reference, coord, layout){

    # set inherited variables
    options <- settings$Plot_Options()
    sourceId <- options$Sample$value
    cellType <- options$Cell_Type$selected
    nPlottedCells <- input$cellsPerPage
    heightPerCell <- options$Cell_Height_Inches$value
    maxCN <- options$Maximum_Copy_Number$value

    # set the working cells
    req(nPlottedCells > 0)
    pageNumber <- as.integer(input$pageNumber)
    plotIndices <- 1:nPlottedCells
    cellIndices <- plotIndices + nPlottedCells * (pageNumber - 1)

    # set axis limits, height and margins
    ylim <- c(0, nPlottedCells)    
    ylim_cell <- c(-0.25, maxCN + 0.25)
    padding <- padding(settings, layout)
    height <- nPlottedCells * heightPerCell + padding$total    
    mai <- NULL

    # set the scCNV data source
    updateSpinnerMessage(session, message = "loading sample data")
    sample <- loadScCnvSample(sourceId)
    setkey(sample$colData, "cell_id")
    # sample$chromEnds <- sample$rowRanges[, max(bin_n, na.rm = TRUE), by = chrom][[2]]
    ploidy <- as.integer(sample$env$PLOIDY)

    # use the mdiTrackImage helper function to create the track image
    image <- mdiTrackImage(layout, height, function(...){

        # set the layout
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        layout(matrix(plotIndices, ncol = 1))

        # plot each cell from the working page
        updateSpinnerMessage(session, message = "plotting cells")
        for(i in cellIndices){

            # initial individual cell plot
            cell_id <- sample$constants[[paste0(cellType, "_cell_ids")]][i]
            par(
                mai = c(0, mai[2], 0, mai[4]),
                cex = 1
            ) # of each subplot
            plot(
                x = NA,
                y = NA,
                typ = "n",
                bty = "n",
                xlim = coord$range,
                ylim = ylim_cell,
                ylab = paste0("#", cell_id, " CN"),
                xaxs = "i",
                yaxs = "i",
                xaxt = "n"
            )
            rect(coord$start, ploidy - 0.5, coord$end, ploidy + 0.5, col = "grey80")            
            abline(h = 0:maxCN, col = "black")
            abline(h = 0:(maxCN - 1) + 0.5, col = "grey50")
            # abline(v = c(0, sample$chromEnds), col = "grey30")

            # collect the cell's data (if any available)
            if(is.na(cell_id) || is.null(cell_id)) next
            cell <- sample[[paste0(cellType, "Cells")]][[cell_id]]
            window_size <- sample$colData[cell_id, window_size]
            if(is.na(window_size)) next

            # get this cell's window x positions
            ww <- paste("w", window_size, sep = "_")
            bin_wr <- sample$rowRanges[, .SD, .SDcols = ww][[1]]
            if(length(bin_wr) == 0) next
            bin_wr <- sample$rowRanges[bin_wr == TRUE, .(chrom, start, end)]
            if(nrow(bin_wr) == 0) next
            plotted <- bin_wr[,
                chrom == coord$chromosome &
                start <= coord$end &
                end   >= coord$start
            ]
            x <- bin_wr[plotted == TRUE, start + (end - start) / 2]

            # plot the copy number and HMM fit
            points(
                x,
                cell$cn[plotted],
                pch = 16,
                cex = settings$Point_Size$value,
                col = "black"           
            )
            if(!is.null(cell$hmm)) lines(
                x,
                cell$hmm[plotted],
                lwd = 2,
                col = "red3"         
            )  
        }
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
click.scCNV_cellsTrack <- function(track, x, y){
    # custom actions
}
hover.scCNV_cellsTrack <- function(track, x, y){
    # custom actions
}
brush.scCNV_cellsTrack <- function(track, x1, y1, x2, y2){
    # custom actions
}

# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.scCNV_cellsTrack <- function(settings, session, input, reference, track){
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
observers.scCNV_cellsTrack <- list()
navigation.scCNV_cellsTrack <- function(settings, session){

    # handle nav observers, e.g., buttons
    for(x in observers.scCNV_cellsTrack) x$destroy()
    input <- session$input
    observers.scCNV_cellsTrack$prevPage <<- observeEvent(input$prevPage, {
        pageNumber <- as.integer(input$pageNumber)
        if(pageNumber == 1) return()
        updateTextInput(session, "pageNumber", value = pageNumber - 1)
    })
    observers.scCNV_cellsTrack$nextPage <<- observeEvent(input$nextPage, {
        pageNumber <- as.integer(input$pageNumber)
        # if(pageNumber == 1) return()
        updateTextInput(session, "pageNumber", value = pageNumber + 1)
    })

    # return the navigation UI elements
    tagList(
        tags$div(
            class = "trackBrowserInput",
            numericInput(session$ns('cellsPerPage'), "Cells Per Page", value = 10, width = "100px"),
        ),
        tags$div(
            class = "trackBrowserInput",
            actionButton(session$ns('prevPage'), "<"),
        ),
        tags$div(
            class = "trackBrowserInput",
            textInput(session$ns('pageNumber'), "Page", value = 1, width = "41px"),
        ),
        tags$div(
            class = "trackBrowserInput",
            actionButton(session$ns('nextPage'), ">"),
        )
    )
}
