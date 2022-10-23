#----------------------------------------------------------------------
# server components for the adjustCells appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
adjustCellsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'adjustCells'   
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
sample <- extractDataReactive(sourceId)
userOverrides <- reactiveValues()
isUserOverride <- Vectorize(function(cell_id) {
    sourceId <- sourceId()
    !is.null(userOverrides[[sourceId]][[cell_id]]) && userOverrides[[sourceId]][[cell_id]]
})
userModalCN <- reactiveValues()
getUserModalCN <- Vectorize(function(cell_id) {
    sourceId <- sourceId()
    x <- userModalCN[[sourceId]][[cell_id]]
    if(is.null(x)) sample()$cells[[cell_id]]$modal_CN else x
})

# calculate the linear slope of the NR vs. GC plot for helping find replicating cells
# slope is used as the sort parameter for cell plot stacks
slopes <- reactive({
    sample <- sample()
    req(sample)    
    sapply(sample$colData$cell_id, function(cell_id){
        cell <- sample$cells[[cell_id]]
        if(is.null(cell$fit)) return(NA)
        gc <- cell$fit$gcFractions
        i <- gc >= 0.35 & gc <= 0.5
        gc <- gc[i]
        mu <- cell$fit$peak[i]
        coef(lm(mu ~ gc))[2]
    })
})
cells <- reactive({
    sample <- sample()
    req(sample)  
    i <- sample$colData[, {
        switch(input$cellType,
            autoKeep   = !rejected & !isUserOverride(cell_id),
            autoReject =  rejected & !isUserOverride(cell_id),
            userKeep   =  rejected &  isUserOverride(cell_id),
            userReject = !rejected &  isUserOverride(cell_id)
        )
    }]
    n <- sum(i)
    list(
        i = i,
        n = n,
        maxPage = ceiling(n / input$cellsPerPage)
    )
})

#----------------------------------------------------------------------
# a stack of individual cell plots
#----------------------------------------------------------------------
getCellCompositePlot <- function(sample, cell_id){
    fileName <- paste("*", cell_id, "qc", "png", sep = ".")
    pngFile <- Sys.glob(file.path(sample$qcPlotsDir, fileName))
    if(file.exists(pngFile)) tags$img(src = pngFileToBase64(pngFile)) else ""
}
getCellSummary <- function(cell){
    prefix <- session$ns("")
    override <- input$cellType == "autoReject" || input$cellType == "autoKeep"
    label <- if(cell$rejected){
        if(input$cellType == "autoReject") "Keep" else "Reject"
    } else {
        if(input$cellType == "autoKeep") "Reject" else "Keep"
    }
    inputId <- paste("cellSetModalCN", cell$cell_id, sep = "-")
    tags$div(
        style = "display: inline-block; padding: 5px; margin-left:5px;",
        tags$div(paste("cell:", cell$cell_id)),
        tags$div(paste("pass:", cell$pass)),
        tags$div(paste0("window_size: ", cell$window_size, ' (', commify(cell$window_size * 20), 'kb )')),
        tags$div(
            style = "margin-top: 5px;",
            tags$button(
                onclick = paste0("cellKeepReject('", prefix, "', '", cell$cell_id, "', '", override, "')"),
                label
            ),
            tags$input(
                id = inputId,
                type = "number",
                min = 1,
                max = 10,
                step = 1,
                value = isolate({ getUserModalCN(cell$cell_id) }),
                onchange = paste0("cellSetModalCN('", prefix, "', '", cell$cell_id, "', '", inputId, "')"),
            )            
        )
    )
}
output$cellPlots <- renderUI({
    slopes <- slopes()
    req(slopes)       
    cells <- cells()
    req(cells)
    if(cells$n == 0) return("no cells to plot")
    startSpinner(session, message = "loading cell plots")
    sample <- sample()
    colData <- sample$colData[cells$i][order(-slopes[cells$i])]
    cellIndices <- 1:input$cellsPerPage + input$cellsPerPage * (as.integer(input$pageNumber) - 1)
    x <- lapply(cellIndices, function(i){
        if(i > cells$n) return(NULL)
        cell_id <- colData[i, cell_id]
        cell <- sample$cells[[cell_id]]
        tags$div(
            style = paste("border: 1px solid grey; white-space: nowrap;"),
            if(cell$pass == 0) tags$div(
                style = "padding: 10px;",
                paste("cell", cell_id, "was rejected prior to normalization due to insufficient read depth")
            ) else tagList(
                getCellCompositePlot(sample, cell_id),
                getCellSummary(as.list(colData[i]))                
            )
        )
    })
    stopSpinner(session)
    x
})

#----------------------------------------------------------------------
# pagination control
#----------------------------------------------------------------------
output$nCellsInType <- renderText({
    cells <- cells()
    req(cells)
    paste0(cells$n, " cells in category (", cells$maxPage, " pages)")
})
observeEvent(input$prevPage, {
    pageNumber <- as.integer(input$pageNumber)
    if(pageNumber == 1) return()
    updateTextInput(session, "pageNumber", value = pageNumber - 1)
})
observeEvent(input$nextPage, {
    pageNumber <- as.integer(input$pageNumber)
    if(pageNumber == cells()$maxPage) return()
    updateTextInput(session, "pageNumber", value = pageNumber + 1)
})

#----------------------------------------------------------------------
# user override control
#----------------------------------------------------------------------
observeEvent(input$cellKeepReject, {
    sourceId <- sourceId()
    userOverrides[[sourceId]][[input$cellKeepReject$cell_id]] <- as.logical(input$cellKeepReject$override)
})
observeEvent(input$cellSetModalCN, {
    sourceId <- sourceId()
    userModalCN[[sourceId]][[input$cellSetModalCN$cell_id]] <- input$cellSetModalCN$value
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    bo <- bm$outcomes
    for(sourceId in names(bo$userOverrides)) userOverrides[[sourceId]] <- bo$userOverrides[[sourceId]]
    for(sourceId in names(bo$userModalCN))   userModalCN[[sourceId]]   <- bo$userModalCN[[sourceId]]
    # settings$replace(bm$settings)
    # updateTextInput(session, 'cellsPerPage', value = bm$input$cellsPerPage)
    # updateTextInput(session, 'pageNumber',   value = bm$input$pageNumber)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(
        userOverrides = userOverrides,
        userModalCN   = userModalCN
    ),
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
