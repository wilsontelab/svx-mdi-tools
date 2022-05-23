#----------------------------------------------------------------------
# reactive components to plot summary results over all selected samples
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
analyzeSNVsServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'analyzeSNVs' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    templates = list(
        file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), 
        file.path(app$sources$suiteGlobalDir, "settings", "svCapture_filters.yml"),
        id
    ),
    fade = FALSE
)
globalSettingsDir <- file.path(app$sources$suiteGlobalDir, "settings")
alignmentSettings <- settingsServer(
    id = 'alignmentSettings',
    parentId = id,
    templates = file.path(globalSettingsDir, "junction_alignment.yml"),
    fade = FALSE,
    title = "Junction Alignment Settings",
    immediate = TRUE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues()

#----------------------------------------------------------------------
# parse filtered SVs and evidence molecules from selected sample(s) and SV(s)
#----------------------------------------------------------------------
# the set of SV junction passing the query filters
filteredSvs <- reactive({ getGenotypedSvs(settings, sampleSelector) })
# the one working SV the user has clicked on
selectedSv <- reactive({ 
    rowI <- svsTable$rows_selected()
    if(is.null(rowI) || rowI == 0) return(NULL)
    svs <- filteredSvs()
    svs[rowI]
})

# ----------------------------------------------------------------------
# susummary counts and SNV rates
# ----------------------------------------------------------------------
bufferedTableServer(
    id = 'aggregatesTable',
    parentId = id,
    parentInput = input,
    selection = 'single',
    tableData = reactive({ tabulateSmallVariants(filteredSvs, settings) })
)

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
locationsPlot <- staticPlotBoxServer(
    'snvLocations', 
    lines = TRUE,
    title = TRUE,
    margins = TRUE,
    legend = TRUE,
    immediate = TRUE,
    template = read_yaml(file.path(app$sources$suiteGlobalDir, "settings", "variant_location_stats.yml")),
    create = function(...){
        plotSnvsByDistance(filteredSvs, settings, locationsPlot)
    }
)

# ----------------------------------------------------------------------
# summary table of all filtered SVs from all selected samples
# ----------------------------------------------------------------------
matchThreshold <- reactive({
    svFilters <- settings$Variant_Options()
    req(svFilters)
    if(svFilters$Allow_Reference_Matches$value) SVX$matchTypes$MISMATCH else SVX$matchTypes$REFERENCE
})
svsTable <- filteredSvsTableServer(id, input, filteredSvs, matchThreshold = matchThreshold)

# ----------------------------------------------------------------------
# alignment of junction bases to reference genome and the two source haplotypes
# ----------------------------------------------------------------------
output$junctionAlignment <- renderText({
    sv <- selectedSv()
    req(sv)
    charPerLine <- alignmentSettings$get("Alignment_Settings", "Bases_Per_Line")
    paste(
        parseSide(sv, 1, charPerLine),
        parseSide(sv, 2, charPerLine),
        sep = "<br><hr>"
    )
})

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
        locationsPlot$settings$replace(bm$outcomes$locationsPlotSettings)
        alignmentSettings$replace(bm$outcomes$alignmentSettings)
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
        locationsPlotSettings  = locationsPlot$settings$all_(),
        alignmentSettings = alignmentSettings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
