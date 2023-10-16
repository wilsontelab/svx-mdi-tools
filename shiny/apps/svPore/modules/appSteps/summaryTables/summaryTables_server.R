#----------------------------------------------------------------------
# server components for the summaryTables appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
summaryTablesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'summaryTables'
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
# get the sample data to plot
#----------------------------------------------------------------------
sourceData <- svPoreSampleSelectorServer("data")
str(sourceData)
#----------------------------------------------------------------------
# individual junction edges
#----------------------------------------------------------------------
# junctionEdges <- function(){
#     samples <- sourceData$samples() 
#     req(samples)   
#     sourceData$sourceId()
    
# }
# sourceSamplesReactive <- reactive({ 

#     sourceSamplesFile <- expandSourceFilePath(sourceId, "sourceSamples.rds")
#     if(!file.exists(sourceSamplesFile)){
#         startSpinner(session, message = "reading source samples")
#         x <- readRDS(getSourceFilePath(sourceId, "segmentsFile")) 
#         x <- x[, .(nSegments = .N), keyby = .(sample)]
#         saveRDS(x, sourceSamplesFile)
#         stopSpinner(session)
#     }
#     readRDS(sourceSamplesFile)
# })
junctionClusters <- junctionClustersReactive(sourceData$sourceId, session)
bufferedTableBoxServer(
    "junctionClusters",
    tableData = junctionClusters, # reactive, or function with no arguments, that returns the table data
    selection = 'none',
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
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
