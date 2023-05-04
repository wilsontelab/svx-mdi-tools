#----------------------------------------------------------------------
# server components for the summaryPlots appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
summaryPlotsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'summaryPlots'
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
sourceId <- dataSourceTableServer("source", selection = "single") 
svNodes <- svNodesReactive(sourceId, session)
filteredNodes <- filteredNodesReactive(svNodes, input)

#----------------------------------------------------------------------
# construct the 2D plot of SV edge metrics for finding artifact classes
#----------------------------------------------------------------------
plotNodes <- plotNodesReactive(filteredNodes, input)
summaryPlot <- summaryPlotServer(plotNodes, input)

#----------------------------------------------------------------------
# construct a plot that profile the alignment of a single selected molecule
#----------------------------------------------------------------------
segments <- segmentsReactive(svNodes, plotNodes, summaryPlot)
moleculePlot <- moleculePlotServer(sourceId, segments, input, session)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
        # outcomes <<- listToReactiveValues(bm$outcomes)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = reactive({ list(
    ) }),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
