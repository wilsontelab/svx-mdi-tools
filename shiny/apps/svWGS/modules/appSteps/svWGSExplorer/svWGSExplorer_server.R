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
    legend = TRUE,
    immediate = TRUE,
    create = function(...){

#         # load targets bed, assumed to be the same for all sample sources
#         assignments <- sampleSelector$selectedAssignments()
#         req(assignments)
#         req(nrow(assignments) > 0)
#         sourceId <- assignments[, Source_ID[1]]


        depth <- getGenomeCoverage(sampleSelector)[[2]]
        maxChunkIndices <- depth[, .(maxChunkIndex = max(chunkIndex)), by = chromIndex]
        maxChunkIndices <- maxChunkIndices[order(chromIndex), c(0, cumsum(maxChunkIndex))]

        plot(
            depth$chunkIndex + maxChunkIndices[depth$chromIndex], 
            depth$coverage, pch = "."
        )
        abline(v = maxChunkIndices)


        # # initialize the data
        # svFilters <- settings$SV_Filters()
        # stepSettings <- settings$Plot_Settings()

        # svs <- filteredSvs()[, .(TARGET_POS_1, TARGET_POS_2)]

        # svs[, ':='(
        #     size = abs(TARGET_POS_2 - TARGET_POS_1 + 1),
        #     center = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1 + 1) / 2
        # )]
        # svPointColors <- svPointColors()

#         # initialize the plot
#         par(mar = c(4.1, 4.1, 0.1, 0.1))
#         xlim <- c(min(targets$paddedStartI), max(targets$paddedEndI))
#         ylim <- c(0, if(any(c("tt", "ta", "aa") %in% targetClasses())){
#             targets[, max(paddedEndI)]
#         } else {
#             targets[, .(x = endI - paddedStartI + 1), by = regionName][, max(x)]
#         })
        
#         plot(
#             NA, NA, typ = "n",
#             xlim = xlim,
#             ylim = ylim,
#             xlab = "SV Center (Mbp)",
#             ylab = "SV Size (bp)",
#             xaxs = "i", 
#             yaxs = "i",
#             xaxt = "n"
#         )

#         # shade and demarcate the capture target regions
#         targets[, {
#             xinc <- ylim[2] / 2
#             polygon(
#                 x = c(startI, startI + xinc, endI + xinc, endI, startI), 
#                 y = c(ylim[1], ylim[2], ylim[2], ylim[1], ylim[1]), 
#                 border = NA, col = "grey90"
#             )
#             polygon(
#                 x = c(startI, startI - xinc, endI - xinc, endI, startI), 
#                 y = c(ylim[1], ylim[2], ylim[2], ylim[1], ylim[1]), 
#                 border = NA, col = "grey90"
#             )
#             mtext(paste(regionName, chrom, sep = ","), side = 1, line = 2.25, at = centerI, cex = 1)
#             paddedSize <- paddedEndI - paddedStartI + 1
#             unit <- paddedSize / 10
#             at <- seq(unit, paddedSize - unit, unit)
#             axis(1, at = paddedStartI + at, labels = round((start - (startI - paddedStartI) + at) / 1e6, 2))
#         }, by = regionName]    
#         for(i in unlist(targets[, .(paddedStartI, startI, endI, paddedEndI)])){
#             abline(-i * 2,  2, col = "grey60")
#             abline( i * 2, -2, col = "grey60")
#         }

#         # plot the SV points on top
#         points(
#             svs[, center], 
#             svs[, size], 
#             pch = 20, 
#             cex = stepSettings$Point_Size$value,
#             col = svPointColors$colors
#         )

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
        # locationsPlot$settings$replace(bm$outcomes$locationsPlotSettings)
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
        # locationsPlotSettings  = locationsPlot$settings$all_(),
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
