#----------------------------------------------------------------------
# reactive components to plot summary results for non-SV artifacts
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
artifactExplorerServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'artifactExplorer' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    # templates = list(id),
    fade = FALSE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues()

# ----------------------------------------------------------------------
# plot the placement of Tx class SVs, including presumed artifacts
# ----------------------------------------------------------------------
distancePlot <- staticPlotBoxServer(
    'distancePlot',
    points    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = FALSE,
    create = function(...){
        makeArtifactDistancePlot(settings, sampleSelector, distancePlot)
    }
)

# ----------------------------------------------------------------------
# plot the microhomology profile of presumed SV artifact classses vs. true SVs
# ----------------------------------------------------------------------
uHomPlot <- staticPlotBoxServer(
    'microhomologyPlot',
    legend    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = FALSE,
    create = function(...){
        makeArtifactMicrohomologyPlot(settings, sampleSelector, uHomPlot)  
    }
)

# ----------------------------------------------------------------------
# define bookmarking actions
# ----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    sampleSet <- bm$input[['sampleSelector-sampleSet']]
    sampleSelector$setSampleSet(sampleSet) 
    if(!is.null(bm$outcomes)) {
        outcomes <<- listToReactiveValues(bm$outcomes)
        sampleSelector$setSelectedSamples(sampleSet, bm$outcomes$samples)
        distancePlot$settings$replace(bm$outcomes$distancePlotSettings)
        uHomPlot$settings$replace(bm$outcomes$uHomDistributionSettings)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    samples  = sampleSelector$selectedSamples,
    outcomes = reactive({ list(
        samples = sampleSelector$selectedSamples(),
        distancePlotSettings = distancePlot$settings$all_(),
        uHomDistributionSettings = uHomPlot$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
