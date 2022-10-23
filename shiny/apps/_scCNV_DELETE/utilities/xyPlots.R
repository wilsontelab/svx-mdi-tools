#----------------------------------------------------------------------
# create a stack of XY scatterplots of individual cells
#----------------------------------------------------------------------
plot_scCNV_cells <- function(input, settings){


    # # set inherited variables
    # settings <- settings$Plot_Options()
    # cellType <- input$cellType
    # nPlottedCells <- input$cellsPerPage
    # heightPerCell <- settings$Track_Height_Pixels$value
    # maxCN <- settings$Maximum_Copy_Number$value
    # width <- settings$Plot_Width_Pixels$value
    # height <- nPlottedCells * heightPerCell

    # # set the working cells
    # pageNumber <- as.integer(input$pageNumber)
    # plotIndices <- 1:nPlottedCells
    # cellIndices <- plotIndices + nPlottedCells * (pageNumber - 1)

    # # set margins
    # axisMarginInches   <- 1
    # nullMarginInches   <- 0.05
    # bottomMarginInches <- axisMarginInches 
    # leftMarginInches   <- axisMarginInches
    # topMarginInches    <- nullMarginInches
    # rightMarginInches  <- nullMarginInches

    # # set axis limits
    # xlim <- range(x$rowRanges$bin_n, na.rm = TRUE)
    # ylim <- c(-0.25, maxCN + 0.25)

    # # set inherited variables
    # settings <- settings$Plot_Options()    
    # setVarValue <- function(key){
    #     if(is.null(input[[key]])) settings[[key]]value else input[[key]]
    # }
    # cellType <- setVarValue("Cell_Type")
    # nPlottedCells <- setVarValue("Cells_Per_Page")


    # pageNumber <- as.integer(input$pageNumber)


    # heightPerCell <- settings$Track_Height_Pixels$value
    # maxCN <- settings$Maximum_Copy_Number$value
    # width <- settings$Plot_Width_Pixels$value
    # height <- nPlottedCells * heightPerCell

    # # set derived values
    # plotIndices <- 1:nPlottedCells
    # cellIndices <- plotIndices + nPlottedCells * (pageNumber - 1)

    # # initialize composite plot
    # layout(matrix(plotIndices, ncol = 1))

    # # plot each cell from the working page
    # for(i in cellIndices){

    #     # initial individual cell plot
    #     par(
    #         mai = c(nullMarginInches, leftMarginInches, nullMarginInches, rightMarginInches), 
    #         cex = 1 / 0.66
    #     ) # of each subplot
    #     plot(
    #         x = NA,
    #         y = NA,
    #         typ = "n",
    #         xlim = xlim,
    #         ylim = ylim,
    #         ylab = "CN",
    #         xaxs = "i",
    #         yaxs = "i",
    #         xaxt = "n",
    #         bty = "n"
    #     )
    #     abline(h = 0:maxCN, col = "black")
    #     abline(h = 0:(maxCN - 1) + 0.5, col = "grey50")
    #     abline(v = c(0, x$chromEnds), col = "grey30")

    #     # collect the cell's data (if any available)
    #     cell_id <- x$constants[[paste0(cellType, "_cell_ids")]][i]
    #     if(is.na(cell_id) || is.null(cell_id)) next
    #     cell <- x[[paste0(cellType, "Cells")]][[cell_id]]
    #     window_size <- x$colData[cell_id, window_size]
    #     if(is.na(window_size)) next

    #     # get this cell's window x positions
    #     ww <- paste("w", window_size, sep = "_")
    #     bin_wr <- x$rowRanges[, .SD, .SDcols = ww][[1]]
    #     bin_n <- x$rowRanges[bin_wr == TRUE, bin_n]

    #     # plot the copy number and HMM fit
    #     points(
    #         bin_n,
    #         cell$cn,
    #         pch = 16,
    #         cex = settings$Point_Size$value,
    #         col = "black"           
    #     )
    #     if(!is.null(cell$hmm)) points(
    #         bin_n,
    #         cell$hmm,
    #         pch = 16,
    #         cex = 0.5,
    #         col = "red3"         
    #     )  
    # }

    # # send out results to mdiInteractivePlot
    # list(
    #     layout = list(
    #         width = width,
    #         height = height,
    #         dpi = 96,
    #         mai = c(bottomMarginInches, leftMarginInches, topMarginInches, rightMarginInches), # of whole plot
    #         xlim = xlim,
    #         ylim = c(0, nPlottedCells)
    #     )
    #     # ,
    #     # parseLayout = function(x, y) list(x, y, layout) # to convert to plot space in a multi-plot layout
    # )
}
