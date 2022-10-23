#----------------------------------------------------------------------
# server components for the markCells appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
markCellsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'markCells'
fileName <- paste0(app$NAME, "-", id, "-", module, ".png")     
pngFile <- file.path(sessionDirectory, fileName)
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# sample selection and loading
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
projectName <- projectNameReactive(sourceId)
sampleData <- sampleDataReactive(sourceId)

#----------------------------------------------------------------------
# individual cell plots
#----------------------------------------------------------------------
cellPlots <- reactive({
    x <- sampleData()
    req(x)
    startSpinner(session, message = "plotting cells")

    # set inherited variables
    settings <- settings$Plot_Options()
    cellType <- input$cellType
    nPlottedCells <- input$cellsPerPage
    heightPerCell <- settings$Track_Height_Pixels$value
    maxCN <- settings$Maximum_Copy_Number$value
    width <- settings$Plot_Width_Pixels$value
    height <- nPlottedCells * heightPerCell

    # set the working cells
    pageNumber <- as.integer(input$pageNumber)
    plotIndices <- 1:nPlottedCells
    cellIndices <- plotIndices + nPlottedCells * (pageNumber - 1)

    # set margins
    axisMarginInches   <- 1
    nullMarginInches   <- 0.05
    bottomMarginInches <- axisMarginInches 
    leftMarginInches   <- axisMarginInches
    topMarginInches    <- nullMarginInches
    rightMarginInches  <- nullMarginInches

    # set axis limits
    xlim <- range(x$rowRanges$bin_n, na.rm = TRUE)
    ylim <- c(-0.25, maxCN + 0.25)

    # initialize composite plot
    png(
        pngFile,
        width  = width,
        height = height,
        units  = "px",
        pointsize = 6,
        res = 96,
        type = "cairo"
    )
    layout(matrix(plotIndices, ncol = 1))

    # plot each cell from the working page
    for(i in cellIndices){

        # initial individual cell plot
        par(mai = c(nullMarginInches, leftMarginInches, nullMarginInches, rightMarginInches), cex = 1 / 0.66) # of each subplot
        plot(
            x = NA,
            y = NA,
            typ = "n",
            xlim = xlim,
            ylim = ylim,
            # xlab = if(i == nPlottedCells) "Genome Bin Index (20 kb bins)" else NULL,
            ylab = "CN",
            xaxs = "i",
            yaxs = "i",
            xaxt = "n",
            bty = "n"
        )
        abline(h = 0:maxCN, col = "black")
        abline(h = 0:(maxCN - 1) + 0.5, col = "grey50")
        abline(v = c(0, x$chromEnds), col = "grey30")

        # collect the cell's data (if any available)
        cell_id <- x$constants[[paste0(cellType, "_cell_ids")]][i]
        if(is.na(cell_id) || is.null(cell_id)) next
        cell <- x[[paste0(cellType, "Cells")]][[cell_id]]
        window_size <- x$colData[cell_id, window_size]
        if(is.na(window_size)) next

        # get this cell's window x positions
        ww <- paste("w", window_size, sep = "_")
        bin_wr <- x$rowRanges[, .SD, .SDcols = ww][[1]]
        bin_n <- x$rowRanges[bin_wr == TRUE, bin_n]

        # plot the copy number and HMM fit
        points(
            bin_n,
            cell$cn,
            pch = 16,
            cex = settings$Point_Size$value,
            col = "black"           
        )
        if(!is.null(cell$hmm)) points(
            bin_n,
            cell$hmm,
            pch = 16,
            cex = 0.5,
            col = "red3"         
        )  
    }
    dev.off()
    stopSpinner(session)

    # send out results to mdiInteractivePlot
    list(
        pngFile = pngFile,
        layout = list(
            width = width,
            height = height,
            dpi = 96,
            mai = c(bottomMarginInches, leftMarginInches, topMarginInches, rightMarginInches), # of whole plot
            xlim = xlim,
            ylim = c(0, nPlottedCells)
        )
        # ,
        # parseLayout = function(x, y) list(x, y, layout) # to convert to plot space in a multi-plot layout
    )
})
mdiInteractivePlotServer(
    "cellPlots",   
    hover = FALSE,    
    click = FALSE,
    brush = FALSE,
    delay = 500,
    contents = cellPlots
)

#----------------------------------------------------------------------
# pagination control
#----------------------------------------------------------------------
observeEvent(input$prevPage, {
    pageNumber <- as.integer(input$pageNumber)
    if(pageNumber == 1) return()
    updateTextInput(session, "pageNumber", value = pageNumber - 1)
})
observeEvent(input$nextPage, {
    pageNumber <- as.integer(input$pageNumber)
    # if(pageNumber == 1) return()
    updateTextInput(session, "pageNumber", value = pageNumber + 1)
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

# Classes ‘data.table’ and 'data.frame':  120949 obs. of  35 variables:
#  $ chrom      : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start      : int  1 20001 40001 60001 80001 100001 120001 140001 160001 180001 ...
#  $ end        : int  20000 40000 60000 80000 100000 120000 140000 160000 180000 200000 ...
#  $ gc_fraction: num  0 0 0 0 0 0 0 0 0 0 ...
#  $ mappability: num  0 0 0 0 0 0 0 0 0 0 ...
#  $ autosome   : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ chrom_bin_n: int  1 2 3 4 5 6 7 8 9 10 ...
#  $ bad_region : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ bin_n      : int  1 2 3 4 5 6 7 8 9 10 ...
#  $ w_1        : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ w_3        : logi  FALSE TRUE FALSE FALSE TRUE FALSE ...
#  - attr(*, ".internal.selfref")=<externalptr>
# List of 5
#  $ 0  :Classes ‘data.table’ and 'data.frame':   7115 obs. of  4 variables:      
#   ..$ cn : num [1:7115] NA NA NA NA NA ...
#   ..$ hmm: int [1:7115] 2 2 2 2 2 2 2 2 2 2 ...
#   ..$ cnc: num [1:7115] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ cnv: int [1:7115] NA NA NA NA NA NA NA NA NA NA ...
#   ..- attr(*, ".internal.selfref")=<externalptr>
