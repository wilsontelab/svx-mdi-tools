# support for artifactExplorer plots

#----------------------------------------------------------------------
# establish override settings common to all artifact molecules (mimicking svx_filters.yml)
#----------------------------------------------------------------------
artifactSettings <- list( 
    SV_Filters = list(
        Min_SV_Size = list(value = 2000), # exclude only the smallest "SVs" from plots
        Max_SV_Size = list(value = 1e9),  # i.e., no limit
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
distancePlotSettings <- function(settings) list(
    SV_Filters = function(){ c(artifactSettings$SV_Filters, list(
        SV_Type = list(value = c("Del", "Dup", "Inv", "Trans")), 
        Min_Split_Reads = list(value = 0) # allow gap molecules in the distance plot
    )) },
    Capture_Filters = function(){ c(artifactSettings$Capture_Filters, list(
        Min_Read_Count = list(value = 3), # thus, distance plot suppresses chimeric PCR, we are after ligation artifacts # nolint
        Duplex_Filter = list(value = settings$Parameters()$Duplex_Filter$value)
    )) }
)
# ----------------------------------------------------------------------
makeArtifactDistancePlot <- function(settings, sampleSelector, plot){
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
        distancePlotSettings(settings), 
        sampleSelector, 
        isCapture = TRUE,
        targetClasses = distancePlotTargetClasses
    ) 

    # calculate SV counts by target class genome wide, then discard translocations for plot
    classCounts <- svs[, .N, by = TARGET_CLASS]
    classCounts <- paste(
        paste("Net Coverage =", commify(round(onTargetCoverage))), 
        paste("#SVs On Target =", classCounts[TARGET_CLASS != "t-", commify(sum(N))]),
        paste("#SVs Off Target =", classCounts[TARGET_CLASS == "t-", commify(sum(N))]),
        sep = "\n"
    )
    svs <- svs[edgeType != svx_edgeTypes$TRANSLOCATION]

    # calculate distance of the farthest SV junction point from the relevant target region center
    svs <- svs[TARGET_REGION %in% targetRegionNames, .( # must be a single target region
        distanceBin = {
            x <- c(POS_1, POS_2) - targetRegionCenters[[TARGET_REGION]]
            x <- x[which.max(abs(x))]
            as.integer(abs(x) / binSize) * binSize * sign(x)
        }
    ), by = c("PROJECT", "SV_ID")] 

    # calculate the binned Tx SV count, normalized to TT Proper
    x <- svs[, 
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
    plot$initializeFrame(
        xlim = xlim, 
        ylim = ylim,
        xlab = "Distance of 2nd End from Capture Target (bp)",
        ylab = "Normalized SV Count"
    )
    abline(v = 0, col = CONSTANTS$plotlyColors$black)       
    isZero <- x$normalizedSvCount == 0       
    abline(h = 0, col = CONSTANTS$plotlyColors$green) 
    abline(h = max(x$normalizedSvCount), col = CONSTANTS$plotlyColors$blue)
    plot$addLines(
        x = x$distanceBin,
        y = x$normalizedSvCount,
        col = CONSTANTS$plotlyColors$grey
    )
    plot$addPoints(
        x = x$distanceBin[!isZero],
        y = x$normalizedSvCount[!isZero],
        col = CONSTANTS$plotlyColors$blue
    )  
    # abline(h = mean(x$normalizedSvCount), col = CONSTANTS$plotlyColors$red) 
    text(
        min(xlim), 
        max(ylim) * 0.95, 
        labels = classCounts,
        adj = c(0, 1)
    )
}

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
uHomCaptureSettings <- function(settings){
    c(artifactSettings$Capture_Filters, list(
        Duplex_Filter = list(value = settings$Parameters()$Duplex_Filter$value)
    ))
}
uHomPlotSettings <- function(settings) list(
    ligationArtifact = list(
        SV_Filters = function(){ 
            c(uHomSvSettings(), list(
                SV_Type = list(value = c("Del", "Dup", "Inv", "Trans"))
            ))
        },
        Capture_Filters = function(){ 
            c(uHomCaptureSettings(settings), list(
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
            c(uHomCaptureSettings(settings), list(
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
            c(uHomCaptureSettings(settings), list(
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
makeArtifactMicrohomologyPlot <- function(settings, sampleSelector, plot){
    param <- settings$Parameters()
    maxMH <- param$Max_Microhomology$value
    uHomPlotSettings <- uHomPlotSettings(settings)
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
    plot$initializeFrame(
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
        plot$addLines(
            x = x$MICROHOM_LEN,
            y = x$freq,
            col = artifactTypeColors[[type]]
        )
    })
    plot$addLegend(
        legend = paste0(artifactTypeLabels, " (", counts, ")"),
        col    = unlist(artifactTypeColors)
    ) 
}
