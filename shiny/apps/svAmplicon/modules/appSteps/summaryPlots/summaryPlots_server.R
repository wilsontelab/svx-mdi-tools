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
junctions <- junctionsReactive(samples, moleculeTypes, selectedAmplicons) # filtered by quality thresholds

#----------------------------------------------------------------------
# construct the junction types table
#----------------------------------------------------------------------
junctionTypes <- junctionTypesReactive(junctions)
junctionTypesTable <- junctionTypesTableServer(id, input, junctionTypes, selection = "multiple")
junctionTypesJunctions <- junctionTypesJunctionsReactive(junctions, junctionTypes, junctionTypesTable) # i.e., filtered by table selection

#----------------------------------------------------------------------
# make junction plots
#----------------------------------------------------------------------
junctionPlotData <- junctionPlotDataReactive(selectedAmplicons, junctionTypesJunctions)
svTrianglePlot <- svTrianglePlotServer(junctionPlotData)
positionDensityPlot <- positionDensityPlotServer(junctionPlotData)
sizeDensityPlot <- sizeDensityPlotServer(junctionPlotData)

#----------------------------------------------------------------------
# construct the filtered junctions table
#----------------------------------------------------------------------
junctionsTableServer(id, input, junctionPlotData)

#----------------------------------------------------------------------
# make junction position density plots
#----------------------------------------------------------------------
# moleculeTypeData <- reactive({
#     I <- moleculeTypesTable$rows_selected()
#     req(I)
#     tableFilteredMoleculeTypes()[I]
# })
# moleculeTypeExpansion <- reactive({
#     moleculeTypeData <- moleculeTypeData()
#     cols <- names(moleculeTypeData)
#     moleculeTypeData <- moleculeTypeData[, .SD, .SDcols = cols[!(cols %in% c("dotplotL","dotplotR"))]]
#     data.table(
#         key_ = names(moleculeTypeData), 
#         value_ = as.character(unlist(moleculeTypeData))
#     )        
# })

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
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
