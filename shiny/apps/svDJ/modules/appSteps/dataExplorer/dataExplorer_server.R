#----------------------------------------------------------------------
# server components for the dataExplorer appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
dataExplorerServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'dataExplorer'
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
moduleEnv <- environment()

#----------------------------------------------------------------------
# load the data source and it requested objects
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
sourceHasChanged <- reactiveVal(0)
observeEvent(sourceId(), {
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "load source data objects")
    for(sourceObject in names(options$sourceObjects)) assign(
        sourceObject, 
        readRDS(getSourceFilePath(sourceId, options$sourceObjects[[sourceObject]])), # TODO: extend beyond RDS files
        envir =  moduleEnv
    )
    sourceHasChanged( sourceHasChanged() + 1 )
    invalidateTable( invalidateTable() + 1 )
    invalidatePlot( invalidatePlot() + 1 )
    stopSpinner(session)
})

#----------------------------------------------------------------------
# show requested data object structures
#----------------------------------------------------------------------
output$sourceObjectStructure <- renderPrint({
    sourceId <- sourceId()
    req(sourceId)
    req(input$sourceObjectChoices)
    tryCatch({
        x <- get(input$sourceObjectChoices, envir =  moduleEnv)
        str(x)
    }, error = function(e){
        print(e)
        NULL
    })
})

#----------------------------------------------------------------------
# initialize the code editors
#----------------------------------------------------------------------
tableDataEditorId <- "tableDataEditor"
tableDataEditorCssId <- session$ns(tableDataEditorId)
tableDataEditorContentsId <- paste(tableDataEditorId, "contents", sep = "-")
plotEditorId <- "plotEditor"
plotEditorCssId <- session$ns(plotEditorId)
plotEditorContentsId <- paste(plotEditorId, "contents", sep = "-")
codeBuffer <- reactiveValues()
codeBuffer[[tableDataEditorContentsId]] <- "# write code here\n# e.g., my_table[...]\n\n\n"
codeBuffer[[plotEditorContentsId]]      <- "# write code here\n# use `dt` for the data.table above\n\n\n"
initializeCodeEditors <- function(){
    isolate({ 
        startSpinner(session, message = "initializing code editors")
        session$sendCustomMessage("initializeAceCodeEditor", tableDataEditorCssId) 
        session$sendCustomMessage("initializeAceCodeEditor", plotEditorCssId) 
        session$sendCustomMessage("setAceCodeContents", list(
            editorId = tableDataEditorCssId,
            code = codeBuffer[[tableDataEditorContentsId]]
        )) 
        session$sendCustomMessage("setAceCodeContents", list(
            editorId = plotEditorCssId,
            code = codeBuffer[[plotEditorContentsId]] 
        )) 
        stopSpinner(session)
    })    
}
initializeCodeEditors()

#----------------------------------------------------------------------
# construct the data.table from user code
#----------------------------------------------------------------------
invalidateTable <- reactiveVal(0)
observeEvent({
    input$updateTableData
    invalidateTable()
}, {
    session$sendCustomMessage("getAceCodeContents", list(editorId = tableDataEditorCssId))
    session$sendCustomMessage("getAceCodeContents", list(editorId = plotEditorCssId))
    setTimeout(function(...) invalidatePlot( invalidatePlot() + 1 ), delay = 500)    
})
tableData <- reactive({
    editor <- input[[tableDataEditorContentsId]]
    codeBuffer[[tableDataEditorContentsId]] <- editor$code
    tryCatch({
        eval(parse(text = editor$code))
    }, error = function(e) {
        safeError(e)
        NULL
    })
})
bufferedTableServer(
    "tableContents",
    id,
    input,
    tableData,
    selection = 'none',
    options = list()
)

#----------------------------------------------------------------------
# construct the plot from data.table and user code
#----------------------------------------------------------------------
invalidatePlot <- reactiveVal(0)
plotFile <- file.path(sessionDirectory, "svDJ_dataExplorer.png")
plotWidth <- 600
plotHeight <- 435
plotRes <- 96
pointSize <- 10
observeEvent(input$updatePlot, {
    session$sendCustomMessage("getAceCodeContents", list(editorId = plotEditorCssId))
})
plotReactive <- reactive({ 
    invalidatePlot()
    editor <- input[[plotEditorContentsId]]
    unlink(plotFile)    
    if(isTruthy(editor$code)){
        codeBuffer[[plotEditorContentsId]] <- editor$code
        dt <- tableData() # plot may or may not depend on the customized data.table
        png(
            filename = plotFile,
            width = plotWidth, height = plotHeight, units = "px", pointsize = pointSize,
            res = plotRes, type = "cairo"
        )
        tryCatch({
            eval(parse(text = editor$code))
        }, error = function(e) print(e))
        dev.off()          
    }
    if(!file.exists(plotFile)) NULL else list(
        pngFile = plotFile,
        layout = list(
            width = plotWidth,
            height = plotHeight,
            pointsize = pointSize,
            dpi = plotRes
        )
    )
})
mdiInteractivePlotServer(
    "plot",   
    hover = FALSE,    
    click = FALSE,
    brush = FALSE,
    delay = 500,
    contents = plotReactive
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    for(key in names(codeBuffer)){
        if(!is.null(bm$outcomes[[key]])) codeBuffer[[key]] <- bm$outcomes[[key]]
    }
    initializeCodeEditors()
    invalidateTable()
    invalidatePlot()
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = codeBuffer,
    isReady = reactive({ getStepReadiness(options$source) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
