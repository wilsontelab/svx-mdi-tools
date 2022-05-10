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
# settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
#     id = 'settings',
#     parentId = id,
#     templates = list(

#             # file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), 
#         # file.path(app$sources$suiteGlobalDir, "settings", "svCapture_filters.yml"),    id
#     ),
#     fade = FALSE
# )
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues()

#----------------------------------------------------------------------
# generate the list of all filtered SVs from all selected samples
#----------------------------------------------------------------------
# filteredSvs <- reactive({ 
#     getFilteredSvs(settings, sampleSelector, isCapture = TRUE) 
# })

# ----------------------------------------------------------------------
# figures related to t- SVs
# ----------------------------------------------------------------------
distancePlot <- staticPlotBoxServer(
    'distancePlot',
    legend    = TRUE,
    points    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){

        # load targets bed, assumed to be the same for all sample sources
        assignments <- sampleSelector$selectedAssignments()
        req(assignments)
        req(nrow(assignments) > 0)
        sourceId <- assignments[, Source_ID[1]]
        targetsBed <- loadPersistentFile(sourceId = sourceId, contentFileType = "targetsBed") 
        targets <- persistentCache[[targetsBed]]$data
        targets[, center := start + size / 2]
        targetRegionNames <- c(targets$regionName, paste0("*,", targets$regionName))
        targetRegionCenters <- as.list(c(targets$center, targets$center))
        names(targetRegionCenters) <- targetRegionNames

        # TODO: expose as settings
        binSize <- 50000
        maxBins <- 100
        maxDistance <- binSize * maxBins

        # get the summed onTargetCoverage, i.e., TT Proper, for the libraries being analyzed, for normalization
        libraryMetrics <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
            manifestFile <- getSourceFilePath(sourceId, "manifestFile")
            x <- fread(manifestFile)
            x[, ':='(uniqueId = paste(Project, Sample_ID, sep = ":"))]
            x
        }))
        onTargetCoverage <- libraryMetrics[uniqueId %in% assignments[, uniqueId], sum(onTargetCoverage)]

        # load the relevant SVs, with override of the the target class filter for this plot
        svs <- getFilteredSvs(
            settings, 
            sampleSelector, 
            isCapture = TRUE, 
            targetClasses = SVX$targetClasses$distancePlot # TT, TA, t- types
        ) 

        # calculate distance of the farthest SV junction point from the relevant target region center
        svs <- svs[
            CHROM_1 == CHROM_2 & # only del/dup/inv are informative here
            TARGET_REGION %in% targetRegionNames, 
        .(
            distanceBin = {
                distances <- c(POS_1, POS_2) - targetRegionCenters[[TARGET_REGION]]
                as.integer(distances[which.max(abs(distances))] / binSize) * binSize
            }
        ), by = SV_ID]

        # calculate the binned Tx SV count, normalized to TT Proper
        x <- svs[
            abs(distanceBin) <= maxDistance, 
            .(normalizedSvCount = .N / onTargetCoverage), 
            by = distanceBin
        ]
        # zeros <- data.table(distanceBin = seq(-maxDistance, maxDistance, binSize))
        x <- merge(
            x, 
            data.table(distanceBin = seq(-maxDistance, maxDistance, binSize)
        ), all.y = TRUE)[order(distanceBin)]
        x[is.na(normalizedSvCount), normalizedSvCount := 0]

        # construct the plot
        distancePlot$initializeFrame(
            xlim = c(-1, 1) * maxDistance, 
            ylim = c(0, max(x$normalizedSvCount)),
            xlab = "Distance of Farthest End from Target (bp)",
            ylab = "Normalized Bin SV Count"
        )
        abline(v = 0, col = CONSTANTS$plotlyColors$black)
        isZero <- x$normalizedSvCount == 0        
        distancePlot$addLines(
            x = x$distanceBin,
            y = x$normalizedSvCount,
            col = CONSTANTS$plotlyColors$grey
        )
        # distancePlot$addPoints(
        #     x = x$distanceBin[isZero],
        #     y = x$normalizedSvCount[isZero],
        #     col = CONSTANTS$plotlyColors$grey
        # )
        distancePlot$addPoints(
            x = x$distanceBin[!isZero],
            y = x$normalizedSvCount[!isZero],
            col = CONSTANTS$plotlyColors$blue
        )  
    }
)
# ligationArtifactPlot <- staticPlotBoxServer(
#     'ligationArtifactPlot',
#     legend    = TRUE,
#     points    = TRUE,
#     lines     = TRUE,
#     title     = TRUE,
#     margins   = TRUE,
#     immediate = TRUE,
#     create = function(...){

#         # load the relevant SVs
#         svs <- getFilteredSvs(
#             settings, 
#             sampleSelector, 
#             isCapture = TRUE, 
#             targetClasses = "t-",
#             noSize = TRUE
#         ) 

#         x <- svs[
#             JXN_BASES != "*",
#             .(freq = .N / nrow(svs)),
#             by = MICROHOM_LEN         
#         ]

#         maxMH <- 20

#         x <- merge(
#             x, 
#             data.table(MICROHOM_LEN = seq(-maxMH, maxMH, 1)
#         ), all.y = TRUE)[order(MICROHOM_LEN)]
#         x[is.na(freq), freq := 0]

#         # construct the plot
#         distancePlot$initializeFrame(
#             xlim = c(-1, 1) * maxMH, 
#             ylim = c(0, max(x$freq)) * 1.1,
#             xlab = "Microhomology Length (bp)",
#             ylab = "Frequency"
#         )
#         abline(v = 0, col = CONSTANTS$plotlyColors$black)
#         # isZero <- x$normalizedSvCount == 0        
#         distancePlot$addLines(
#             x = x$MICROHOM_LEN,
#             y = x$freq,
#             col = CONSTANTS$plotlyColors$blue
#         )
#         # distancePlot$addPoints(
#         #     x = x$distanceBin[isZero],
#         #     y = x$normalizedSvCount[isZero],
#         #     col = CONSTANTS$plotlyColors$grey
#         # )
#         # distancePlot$addPoints(
#         #     x = x$distanceBin[!isZero],
#         #     y = x$normalizedSvCount[!isZero],
#         #     col = CONSTANTS$plotlyColors$blue
#         # )  
#     }
# )

# ----------------------------------------------------------------------
# define bookmarking actions
# ----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    sampleSet <- bm$input[['sampleSelector-sampleSet']]
    sampleSelector$setSampleSet(sampleSet) 
    if(!is.null(bm$outcomes)) {
        outcomes <<- listToReactiveValues(bm$outcomes)
        sampleSelector$setSelectedSamples(sampleSet, bm$outcomes$samples)
        # svRatesPlot$settings$replace(bm$outcomes$svRatesPlotSettings)
        # sizeDistribution$settings$replace(bm$outcomes$sizeDistributionSettings)
        # microhomologyDistribution$settings$replace(bm$outcomes$microhomologyDistributionSettings)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    samples  = sampleSelector$selectedSamples,
    outcomes = reactive({ list(
        samples = sampleSelector$selectedSamples()
        # ,
        # svRatesPlotSettings = svRatesPlot$settings$all_(),
        # sizeDistributionSettings = sizeDistribution$settings$all_(),
        # microhomologyDistributionSettings = microhomologyDistribution$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
