#----------------------------------------------------------------------
# reactive components to filter and examine SV locations and junction sequences
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svWGSExplorerServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'svWGSExplorer' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
moduleOptions <- stepModuleInfo[[ app$config$appSteps[[id]]$module ]]
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    templates = list(
        file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), 
        id
    ),
    fade = FALSE,
    size = "m"
)
sampleSelector <- sampleSelectorServer(id = 'sampleSelector', parentId = id)
mapSettings <- mapSettingsServer(id, moduleOptions)
alignmentSettings <- alignmentSettingsServer(id, moduleOptions)
outcomes <- reactiveValues() # logical failure vectors keyed as [[sampleSet]]

#----------------------------------------------------------------------
# parse filtered SVs and evidence molecules from selected sample(s) and SV(s)
#----------------------------------------------------------------------
# the set of SV junction passing the query filters
filteredSvs <- reactive({ getFilteredSvs(settings, sampleSelector) })
# the one working SV the user has clicked on
selectedSv <- reactive({ 
    rowI <- svsTable$rows_selected()
    req(rowI)
    svs <- filteredSvs()
    svs[rowI]
})
# complete supporting molecule evidence on selectedSv()
svMols <- reactive({ getSVMolecules(sampleSelector, selectedSv) })
# map of all base positions at and around the selected SV junction
junctionMap <- reactive({ getJunctionMap(
    svMols(), 
    mapSettings$get("Map_Settings", "Clip_Mode")) 
})
svPointColors <- getSvPointColors(filteredSvs, settings, isCapture = FALSE)

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
locationsPlot <- staticPlotBoxServer(
    'svLocations', 
    points = TRUE,
    legend = TRUE,
    immediate = TRUE,
    create = function(...){
#         # load targets bed, assumed to be the same for all sample sources
#         assignments <- sampleSelector$selectedAssignments()
#         req(assignments)
#         req(nrow(assignments) > 0)
#         sourceId <- assignments[, Source_ID[1]]

        x <- getGenomeCoverage(sampleSelector)

        x[, i := .I] # record I so that gaps will show on screen
        x <- x[x$genmap >= 0.5 & x$excluded <= 65536 / 2]

        normalizedDepth <- x$NA12878 / x$genmap

        modalCN <- 2

        locationsPlot$initializeFrame(
            xlab = "bin",
            ylab = "log2(normalizedDepth)",
            xlim = range(x$i),
            ylim = log2(c(0.75, modalCN + 4) / modalCN)
        )

        verticals <- sapply(unique(x$chrom), function(chr) x[chrom == chr, max(i, na.rm = TRUE)])
        abline(v = c(0, verticals), col = CONSTANTS$plotlyColors$grey)
        abline(h = log2((modalCN + -modalCN:4) / modalCN), col = CONSTANTS$plotlyColors$grey)
        abline(h = 0, col = CONSTANTS$plotlyColors$black)

        locationsPlot$addPoints(
            x = x$i, 
            y = log2(normalizedDepth / median(normalizedDepth))
        )

#         # add a legend
#         pointColorLegend(stepSettings, locationsPlot$settings, svPointColors)
    }
)

#----------------------------------------------------------------------
# SV junction properties plot
#----------------------------------------------------------------------
propertiesPlot <- svPropertiesPlotServer(settings, svPointColors, filteredSvs)

# ----------------------------------------------------------------------
# summary table of all filtered SVs from all selected samples
# ----------------------------------------------------------------------
svsTable <- filteredSvsTableServer(id, input, filteredSvs)

# ----------------------------------------------------------------------
# expanded views of the selected junction
# ----------------------------------------------------------------------
junctionMapServer(output, junctionMap, mapSettings)
junctionNodesPlotServer(svMols)
junctionAlignmentServer(output, junctionMap, alignmentSettings)

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
        propertiesPlot$settings$replace(bm$outcomes$propertiesPlotSettings)
        mapSettings$replace(bm$outcomes$mapSettings)
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
        propertiesPlotSettings = propertiesPlot$settings$all_(),
        mapSettings = mapSettings$all_(),
        alignmentSettings = alignmentSettings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
