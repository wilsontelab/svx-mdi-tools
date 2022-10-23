#----------------------------------------------------------------------
# server components for the adjustFit appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
adjustFitServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'adjustFit'
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
    # settings = id, # for step-level settings
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
data <- extractDataReactive(sourceId)
userOverride <- list()
observeEvent(data(), {
    x <- data()
    req(x)     
    userOverride <<- lapply(x$colData$cell_id, function(cell_id) FALSE)
})

# calculate the linear slope of the NR vs. GC plot for helping find replicating cells
coeff <- reactive({
    x <- data()
    req(x)    
    d <- as.data.table(t(sapply(x$colData$cell_id, function(cell_id){
        cell <- x$cells[[cell_id]]
        if(is.null(cell$fit)) return(c(NA, NA))
        gc <- cell$fit$gcFractions
        i <- gc >= 0.35 & gc <= 0.5
        gc <- gc[i]
        mu <- cell$fit$peak[i]
        coef(lm(mu ~ gc))
    })))  
    d[, cell_id := x$colData$cell_id]
    setkey(d, "cell_id")
    d
})

#----------------------------------------------------------------------
# individual cell plots
#----------------------------------------------------------------------
plotByWindow <- function(d, ylab, ylim = NULL){
    plot(1:length(d), d, bty = "n",
        pch = 19, cex = 0.4, 
            col = rgb(0, 0, 0, 0.2),
        xaxt = "n", xlab = NULL, ylab = ylab, ylim = ylim)
}
getCellCompositePlot <- function(cell){

    # set the composite plot layout
    png(
        filename = pngFile,
        width  = 1.5 * 6, 
        height = 1.35, 
        units = "in", 
        pointsize = 7,
        bg = "white",
        res = 96, # i.e., optimized for on-screen display
        type = "cairo"
    )
    layout(matrix(c(c(1,1), rep(c(2,3), 5)), nrow = 2, ncol = 6))

    # plot NR_map_w vs. gw_w
    nb   <- cell$fit
    col  <- if(cell$rejected) "red3" else "blue"
    peak <- nb$peak  
    mu   <- nb$mu
    par(mar= c(4.1, 4.1, 0.1, 1.1), cex = 1)
    plot(cell$gc_w, cell$NR_map_w, xlim = c(0.35, 0.55), 
            pch = 19, cex = 0.4, 
            col = rgb(0, 0, 0, 0.1),
            xlab = "Fraction GC", ylab = "# of Reads")
    lines(
        nb$gcFractions,
        ifelse(peak < 10, mu, peak) * cell$modal_CN,
        lty = 1, lwd = 1.5, col = col
    )
    for(percentile in c(0.05, 0.95)) lines(
        nb$gcFractions, 
        qnbinom(percentile, size = nb$theta, mu = mu) * cell$modal_CN, 
        lty = 3, lwd = 1.5, col = col
    )

    # plot NR_map_w vs. window index, i.e., pre-normalization
    par(mar= c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plotByWindow(cell$NR_map_w, "# Reads")

    # plot CN vs. window index, i.e., post-normalization
    plotByWindow(cell$cn, "CN", ylim = c(0, 4))
    abline(h=0:4, col = "grey")
    if(!is.null(cell$hmm)) lines(1:length(cell$hmm), cell$hmm, col = "red")

    # return an image tag
    dev.off()
    tags$img(src = pngFileToBase64(pngFile))
}
getCellSummary <- function(cell){

    tags$div(
        style = "display: inline-block; padding: 5px; margin-left:5px;",
        tags$div(paste("cell:", cell$cell_id)),
        tags$div(paste("pass:", cell$pass)),
        tags$div(paste("window_size:", cell$window_size)),
        # tags$div(paste("intercept:", coeff()[cell$cell_id, cell_id])),
        # tags$p(paste("gc:", coeff()[cell$cell_id, gc]))
        # # ,
        # tags$div(paste("gc2:", coeff()[cell$cell_id, 3]))
        ""
    )

}
output$cellPlots <- renderUI({
    x <- data()
    req(x)
    # startSpinner(session, message = "plotting cells")

    # set inherited variables
    # settings <- settings$Plot_Options()
    cellType <- input$cellType
    nPlottedCells <- input$cellsPerPage
    # heightPerCell <- settings$Track_Height_Pixels$value
    # maxCN <- settings$Maximum_Copy_Number$value
    # width <- settings$Plot_Width_Pixels$value
    # height <- nPlottedCells * heightPerCell

    pageNumber <- as.integer(input$pageNumber)
    plotIndices <- 1:nPlottedCells
    cellIndices <- plotIndices + nPlottedCells * (pageNumber - 1)

    coeff <- coeff()
    req(coeff)
    userOverride <- unlist(userOverride)
    cellI <- x$colData[, {
        switch(cellType,
            autoKeep   = !rejected & !userOverride,
            autoReject =  rejected & !userOverride,
            userKeep   =  rejected &  userOverride,
            userReject = !rejected &  userOverride
        )
    }]
    colData <- x$colData[cellI][order(-coeff[cellI, gc])]
    lapply(cellIndices, function(i){
        cell_id <- colData[i, cell_id]
        cell <- x$cells[[cell_id]]
        col <- if(cell$rejected) "#cc0000" else "#0000ee"
        tags$div(
            onclick = paste0("setAdjustFitCell('", session$ns(""), "', '", cell_id, "')"),
            style = paste("border: 1px solid", col, "; white-space: nowrap; display: inline-block;"),
            getCellCompositePlot(cell),
            getCellSummary(as.list(colData[i]))
        )
    })
})

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
