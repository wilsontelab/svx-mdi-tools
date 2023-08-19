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
# get the samples and amplicons to plot
#----------------------------------------------------------------------
sampleSet <- sampleSetServer("sampleSet", id)
samples <- samplesReactive(sampleSet)
amplicons <- ampliconsReactive(samples)

#----------------------------------------------------------------------
# construct the amplicon table and cascading reactives
#----------------------------------------------------------------------
ampliconsTable <- ampliconsTableServer(id, input, amplicons, selection = "multiple")
selectedAmplicons <- selectedAmpliconsReactive(amplicons, ampliconsTable)
moleculeTypes <- moleculeTypesReactive(samples, selectedAmplicons)

#----------------------------------------------------------------------
# construct the path class table - only these path classes are used to construct the junctions table
#----------------------------------------------------------------------
pathClasses <- pathClassesReactive(moleculeTypes)
pathClassesTable <- pathClassesTableServer(id, input, pathClasses, selection = "multiple")
pathClassMoleculeTypes <- pathClassMoleculeTypesReactive(moleculeTypes, pathClasses, pathClassesTable, selectedAmplicons)

#----------------------------------------------------------------------
# construct the junction types table
#----------------------------------------------------------------------
junctions <- junctionsReactive(samples, selectedAmplicons, pathClassMoleculeTypes)
junctionTypes <- junctionTypesReactive(junctions)
junctionTypesTable <- junctionTypesTableServer(id, input, junctionTypes, selection = "multiple")
junctionTypesJunctions <- junctionTypesJunctionsReactive(junctions, junctionTypes, junctionTypesTable) # i.e., filtered by table selection

#----------------------------------------------------------------------
# make junction plots
#----------------------------------------------------------------------
junctionPlotData <- junctionPlotDataReactive(selectedAmplicons, junctionTypesJunctions)
# svTrianglePlot <- svTrianglePlotServer(junctionPlotData)
svCirclesPlot <- svCirclesPlotServer(junctionPlotData)
positionDensityPlot <- positionDensityPlotServer(junctionPlotData)
sizeDensityPlot <- sizeDensityPlotServer(junctionPlotData)

#----------------------------------------------------------------------
# construct the filtered junctions table
#----------------------------------------------------------------------
junctionsTableServer(id, input, junctionPlotData)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
        svCirclesPlot$settings$replace(bm$outcomes$svCirclesPlotSettings)
        # outcomes <<- listToReactiveValues(bm$outcomes)
    }
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = reactive({ list(
        svCirclesPlotSettings = svCirclesPlot$settings$all_()
    ) }),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
