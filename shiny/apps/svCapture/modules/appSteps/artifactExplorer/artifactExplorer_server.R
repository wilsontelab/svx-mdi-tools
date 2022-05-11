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

#----------------------------------------------------------------------
# establish override settings common to all artifact molecules (mimicking svx_filters.yml)
#----------------------------------------------------------------------
artifactSettings <- list( 
    SV_Filters = list(
        Min_SV_Size = list(value = 10000), # exclude only the smallest "SVs" from plots
        Max_SV_Size = list(value = 1e9),
        Min_Samples_With_SV  = list(value = 1), # artifacts are always single-source molecules
        Max_Samples_With_SV  = list(value = 1),
        Min_Source_Molecules = list(value = 1),
        Max_Source_Molecules = list(value = 1),
        Include_Clips_In_Total = list(value = FALSE),
        Min_Map_Quality = list(value = 50) # only interested in true alignments that could not be filtered by quality # nolint
    ),
    Capture_Filters = list(
        Max_Frac_Shared_Proper = list(value = 2) # i.e., not used
    )    
)

# ----------------------------------------------------------------------
# plot the placement of Tx class SVs, including presumed artifacts
# ----------------------------------------------------------------------
distancePlotTargetClasses <- c("TT", "TA", "t-")
distancePlotSettings <- list(
    SV_Filters = function(){ c(artifactSettings$SV_Filters, list(
        SV_Type = list(value = c("Del", "Dup", "Inv")), # with above, same-chrom with one end in target
        Min_Split_Reads = list(value = 0) # allow gap molecules in the distance plot
    )) },
    Capture_Filters = function(){ c(artifactSettings$Capture_Filters, list(
        Min_Read_Count = list(value = 3), # thus, distance plot suppresses chimeric PCR, we are after ligation artifacts # nolint
        Duplex_Filter = list(value = settings$Parameters()$Duplex_Filter$value)
    )) }
)
distancePlot <- staticPlotBoxServer(
    'distancePlot',
    points    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){
        param <- settings$Parameters()
        binSize <- param$Genome_Bin_Size$value
        maxBins <- param$Number_of_Bins_from_Center$value
        maxDistance <- binSize * maxBins
        
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
        targetX <- c(-1, 1) * as.integer(max(targets$size) / 2)

        # get the summed onTargetCoverage, i.e., TT Proper, for the libraries being analyzed, for normalization
        libraryMetrics <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
            manifestFile <- getSourceFilePath(sourceId, "manifestFile")
            x <- fread(manifestFile)
            x[, ':='(uniqueId = paste(Project, Sample_ID, sep = ":"))]
            x
        }))
        onTargetCoverage <- libraryMetrics[uniqueId %in% assignments[, uniqueId], sum(onTargetCoverage)]
        reportProgress(paste("onTargetCoverage =", onTargetCoverage))

        # load the relevant SVs
        svs <- getFilteredSvs(
            distancePlotSettings, 
            sampleSelector, 
            isCapture = TRUE,
            targetClasses = distancePlotTargetClasses
        ) 

        # calculate distance of the farthest SV junction point from the relevant target region center
        svs <- svs[TARGET_REGION %in% targetRegionNames, .( # must be a single target region
            distanceBin = {
                x <- c(POS_1, POS_2) - targetRegionCenters[[TARGET_REGION]]
                x <- x[which.max(abs(x))]
                as.integer(abs(x) / binSize) * binSize * sign(x)
            }
        ), by = SV_ID]

        # calculate the binned Tx SV count, normalized to TT Proper
        x <- svs[
            abs(distanceBin) <= maxDistance, 
            .(normalizedSvCount = .N / onTargetCoverage), 
            by = distanceBin
        ]
        x <- merge(
            x, 
            data.table(distanceBin = seq(-maxDistance, maxDistance, binSize)
        ), all.y = TRUE)[order(distanceBin)]
        x[is.na(normalizedSvCount), normalizedSvCount := 0]

        # construct the plot
        xlim <- c(-1, 1) * maxDistance
        ylim <- c(0, max(x$normalizedSvCount))
        distancePlot$initializeFrame(
            xlim = xlim, 
            ylim = ylim,
            xlab = "Distance of 2nd End from Capture Target (bp)",
            ylab = "Normalized Bin SV Count"
        )
        # rect(targetX[1], ylim[1], targetX[2], ylim[2], col = CONSTANTS$plotlyColors$grey, border = NA)
        abline(v = 0, col = CONSTANTS$plotlyColors$black)       
        isZero <- x$normalizedSvCount == 0        
        distancePlot$addLines(
            x = x$distanceBin,
            y = x$normalizedSvCount,
            col = CONSTANTS$plotlyColors$grey
        )
        distancePlot$addPoints(
            x = x$distanceBin[!isZero],
            y = x$normalizedSvCount[!isZero],
            col = CONSTANTS$plotlyColors$blue
        )  
    }
)

# ----------------------------------------------------------------------
# plot the microhomology profile of presumed SV artifact classses vs. true SVs
# ----------------------------------------------------------------------
uHomTargetClasses <- list(
    ligationArtifact = c("t-"), # junctions likely to have arisen pre-library and then captured
    chimericArtifact = c("tt"), # junctions likely to have arisen late by inter-region chimeric PCR
    realSvs = c("TT", "TA")     # junctions at least possible to be real
)
uHomSvSettings <- function(){
    c(artifactSettings$SV_Filters, list(
        Min_Split_Reads = list(value = 1) # require a sequenced junction for microhomology
    ))
}
uHomCaptureSettings <- function(){
    c(artifactSettings$Capture_Filters, list(
        Duplex_Filter = list(value = settings$Parameters()$Duplex_Filter$value)
    ))
}
uHomPlotSettings <- list(
    ligationArtifact = list(
        SV_Filters = function(){ 
            c(uHomSvSettings(), list(
                SV_Type = list(value = c("Del", "Dup", "Inv", "Trans"))
            ))
        },
        Capture_Filters = function(){ 
            c(uHomCaptureSettings(), list(
                Min_Read_Count = list(value = 3) # suppress chimeric PCR from ligation artifacts
            ))
        }
    ),
    chimericArtifact = list(
        SV_Filters = function(){ 
            c(uHomSvSettings(), list(
                SV_Type = list(value = c("Del", "Dup", "Inv", "Trans"))
            ))
        },
        Capture_Filters = function(){ 
            c(uHomCaptureSettings(), list(
                Min_Read_Count = list(value = 1) # suppress ligation artifacts from chimeric PCR
            ))
        }
    ),
    realSvs = list(
        SV_Filters = function(){ 
            delOnly <- settings$Parameters()$True_SV_Types$value == "deletionOnly"
            c(uHomSvSettings(), list(
                SV_Type = list(value = if(delOnly) 
                                       c("Del") else 
                                       c("Del", "Dup", "Inv"))
            ))
        },
        Capture_Filters = function(){ 
            c(uHomCaptureSettings(), list(
                Min_Read_Count = list(value = 3) # suppress chimeric PCR
            ))
        }
    )
)
artifactTypeLabels <- list(
    ligationArtifact = "Ligation",
    chimericArtifact = "Chimeric PCR",
    realSvs          = "True SVs"
)
artifactTypeColors <- list(
    ligationArtifact = CONSTANTS$plotlyColors$blue,
    chimericArtifact = CONSTANTS$plotlyColors$red,
    realSvs          = CONSTANTS$plotlyColors$green
)
uHomPlot <- staticPlotBoxServer(
    'microhomologyPlot',
    legend    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){
        param <- settings$Parameters()
        maxMH <- param$Max_Microhomology$value
        types <- names(uHomPlotSettings)

        # calculate the different microhomology distributions
        dist <- lapply(types, function(type){
            svs <- getFilteredSvs(
                uHomPlotSettings[[type]], 
                sampleSelector, 
                isCapture = TRUE,
                targetClasses = uHomTargetClasses[[type]]
            ) 
            x <- svs[
                JXN_BASES != "*",
                .(
                    N = .N,
                    freq = .N / nrow(svs)
                ),
                by = MICROHOM_LEN         
            ]
            x <- merge(
                x, 
                data.table(MICROHOM_LEN = seq(-maxMH, maxMH, 1)
            ), all.y = TRUE)[order(MICROHOM_LEN)]
            x[is.na(freq), ":="(
                N = 0,
                freq = 0
            )] 
            x[, cumfreq := cumsum(freq)]   
            x[MICROHOM_LEN %between% c(-maxMH, maxMH)]          
        })
        names(dist) <- types
        maxFreq <- max(sapply(types, function(type) max(dist[[type]]$freq)))
        counts <- sapply(types, function(type) sum(dist[[type]]$N))

        # construct the plot
        xlim <- c(-1, 1) * maxMH
        uHomPlot$initializeFrame(
            xlim = xlim, 
            ylim = c(0, maxFreq * 1.05),
            xlab = "Microhomology Length (bp)",
            ylab = "Frequency"
        )
        mtext("insertion", side = 1, line = 2, at = xlim[1], adj = 0, cex = 0.95)
        mtext(expression(paste(mu, "homology")), side = 1, line = 2, at = xlim[2], adj = 1, cex = 0.95) 
        abline(v = 0, col = CONSTANTS$plotlyColors$black)  
        lapply(types, function(type){
            x <- dist[[type]]
            uHomPlot$addLines(
                x = x$MICROHOM_LEN,
                y = x$freq,
                col = artifactTypeColors[[type]]
            )
        })
        uHomPlot$addLegend(
            legend = paste0(artifactTypeLabels, " (", counts, ")"),
            col    = unlist(artifactTypeColors)
        ) 
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
