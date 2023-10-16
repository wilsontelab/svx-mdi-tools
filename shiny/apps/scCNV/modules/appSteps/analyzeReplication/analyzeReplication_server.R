#----------------------------------------------------------------------
# server components for the analyzeReplication appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
analyzeReplicationServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'analyzeReplication'
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
# load sample/cell source data
#----------------------------------------------------------------------
sourceIds <- dataSourceTableServer(
    "source", 
    selection = "multiple"
)
cellsByFractionS <- reactive({
    sourceIds <- sourceIds()
    req(sourceIds) 
    app$keepReject$getReplicatingCells(sourceIds)  
})

#----------------------------------------------------------------------
# stacked FAR plots for all replicating cells
#----------------------------------------------------------------------
zoomChrom <- reactiveVal("chr1")
output$chromPlots <- renderUI({
    cellsByFractionS <- cellsByFractionS()
    req(cellsByFractionS)
    zoomChrom <- zoomChrom()
    req(zoomChrom)
    # invalidateZoomPlots()
    startSpinner(session, message = "loading chrom plots")
    labelRow <- createZoomLabelRow(session, zoomChrom, close = FALSE, commit = FALSE)
    cells <- lapply(1:nrow(cellsByFractionS), function(j){
        sourceId <- cellsByFractionS[j, sourceId]
        cell_id  <- cellsByFractionS[j, cell_id]
        project <- getScCnvProjectData(sourceId)
        cell <- project$cells[[cell_id]]
        plotOneCellUI_FAR(sourceId, project, cell, settings, zoomChrom)        
    })
    isolate({ initGenomePlotClicks(initGenomePlotClicks() + 1) })  
    stopSpinner(session)
    tagList(labelRow, cells)
})
initGenomePlotClicks <- reactiveVal(0)
observeEvent(initGenomePlotClicks(), {
    session$sendCustomMessage("cellPlotsWrapperInit", list(
        prefix = session$ns(""),
        divId = "chromPlotsWrapper"
    ))
    session$sendCustomMessage("cellPlotsWrapperUpdate", list(
        prefix = session$ns(""),
        divId = "chromPlotsWrapper",
        cellsPerPage = nrow(cellsByFractionS()),
        short = TRUE
    ))
}, ignoreInit = TRUE)

#----------------------------------------------------------------------
# handle actions to navigate through chromosomes
#----------------------------------------------------------------------
prevNextZoomChrom <- function(inc){
    project <- getScCnvProjectData(sourceIds()[1])
    chroms <- project$windows[[1]][, unique(chrom)]
    i <- which(chroms == zoomChrom()) + inc
    if(length(i) == 0) i <- 1
    if(i < 1 || i > length(chroms)) return(NULL)
    zoomChrom(chroms[i])
}
observeEvent(input$prevZoomChrom, { prevNextZoomChrom(-1) })
observeEvent(input$nextZoomChrom, { prevNextZoomChrom( 1) })

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    bookmarkObserver$destroy()
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
