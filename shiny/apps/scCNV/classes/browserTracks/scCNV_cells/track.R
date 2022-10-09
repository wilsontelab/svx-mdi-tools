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
    heightPerCell <- options$Cell_Height_Pixels$value
    maxCN <- options$Maximum_Copy_Number$value

    # set the working cells
    req(nPlottedCells > 0)
    pageNumber <- as.integer(input$pageNumber)
    plotIndices <- 1:nPlottedCells
    cellIndices <- plotIndices + nPlottedCells * (pageNumber - 1)

    # set margins
    axisMarginInches   <- 1
    nullMarginInches   <- 0
    bottomMarginInches <- axisMarginInches 
    leftMarginInches   <- axisMarginInches
    topMarginInches    <- nullMarginInches
    rightMarginInches  <- nullMarginInches

    # set axis limits, height and margins
    ylim <- c(0, nPlottedCells)    
    ylim_cell <- c(-0.25, maxCN + 0.25)
    padding <- padding(settings, layout)
    height <- nPlottedCells * heightPerCell /  layout$dpi + padding$total    
    mai <- NULL

    # set the scCNV data source
    # TODO: use function to cache this
    dataFilePath <- getSourceFilePath(sourceId, "normalizeFile")
    sample <- readRDS(dataFilePath)
    setkey(sample$colData, "cell_id")
    sample$chromEnds <- sample$rowRanges[, max(bin_n, na.rm = TRUE), by = chrom][[2]]

    # use the mdiTrackImage helper function to create the track image
    image <- mdiTrackImage(layout, height, function(...){

        # set the layout
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        layout(matrix(plotIndices, ncol = 1))

        # plot each cell from the working page
        for(i in cellIndices){

            # initial individual cell plot
            par(
                mai = c(nullMarginInches, mai[2], nullMarginInches, mai[4]),
                cex = 1
            ) # of each subplot
            plot(
                x = NA,
                y = NA,
                typ = "n",
                bty = "n",
                xlim = coord$range,
                ylim = ylim_cell,
                ylab = "CN",
                xaxs = "i",
                yaxs = "i",
                xaxt = "n"
            )
            abline(h = 0:maxCN, col = "black")
            abline(h = 0:(maxCN - 1) + 0.5, col = "grey50")
            # abline(v = c(0, sample$chromEnds), col = "grey30")

            # collect the cell's data (if any available)
            cell_id <- sample$constants[[paste0(cellType, "_cell_ids")]][i]
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
                pch = 16,
                cex = 0.5,
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

navigation.scCNV_cellsTrack <- function(settings, session){
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
