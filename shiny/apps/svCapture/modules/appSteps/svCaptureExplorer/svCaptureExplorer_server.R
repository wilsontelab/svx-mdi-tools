#----------------------------------------------------------------------
# reactive components to filter and examine SV locations and junction sequences
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svCaptureExplorerServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'svCaptureExplorer' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
moduleOptions <- stepModuleInfo[[ app$config$appSteps[[id]]$module ]]
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'stepSettings',
    parentId = id,
    templates = list(
        file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), 
        file.path(app$sources$suiteGlobalDir, "settings", "svCapture_filters.yml"),
        id,
        file.path(app$sources$suiteGlobalDir, "settings", "external_svs.yml")
    ),
    fade = FALSE,
    size = "m"
)
sampleSelector <- sampleSelectorServer(id = 'sampleSelector', parentId = id)
mapSettings <- mapSettingsServer(id, moduleOptions)
alignmentSettings <- alignmentSettingsServer(id, moduleOptions)
targetClasses <- reactive({ 
    SVX$targetClasses[[settings$Capture_Filters()$Target_Class$value]] 
})
outcomes <- reactiveValues() # logical failure vectors keyed as [[sampleSet]]

#----------------------------------------------------------------------
# parse filtered SVs and evidence molecules from selected sample(s) and SV(s)
#----------------------------------------------------------------------
# the set of SV junction passing the query filters
filteredSvs <- reactive({ getFilteredSvs(settings, sampleSelector, isCapture = TRUE) })
# the one working SV the user has clicked on
selectedSv <- reactive({ 
    rowI <- svsTable$rows_selected()
    if(is.null(rowI) || rowI == 0) return(NULL) # req(rowI)
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
svPointColors <- getSvPointColors(filteredSvs, settings, sampleSelector, isCapture = FALSE)

#----------------------------------------------------------------------
# load a set of external comparator SVs
#----------------------------------------------------------------------
externalSvs <- reactive({
    bedFile <- settings$External_SVs()$SVs_Bed_File$value
    if(is.null(bedFile) || bedFile == "") return(NULL) # not req(...)
    bedFile <- file.path(serverEnv$UPLOADS_DIR, bedFile)
    if(!file.exists(bedFile)) return(NULL)
    x <- fread(bedFile)[, 1:3]
    names(x) <- c("chrom", "start", "end")
    x[, ":="(
        size = end - start,
        center = start + (end - start) / 2
    )]
    x
})

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
locationsPlot <- staticPlotBoxServer(
    'svLocations', 
    title = TRUE,
    margins = TRUE,
    legend = TRUE,
    immediate = TRUE,
    create = function(...){

        # load targets bed, assumed to be the same for all sample sources
        assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
        req(assignments)
        req(nrow(assignments) > 0)
        sourceId <- assignments[, Source_ID[1]]
        targetsBed <- loadPersistentFile(sourceId = sourceId, contentFileType = "targetsBed") 
        targets <- persistentCache[[targetsBed]]$data

        # initialize the data
        filters <- settings$SV_Filters()
        stepSettings <- settings$Plot_Settings()
        svs <- filteredSvs()[, .(TARGET_POS_1, TARGET_POS_2)]
        svs[, ':='(
            size = abs(TARGET_POS_2 - TARGET_POS_1 + 1),
            centerI = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1 + 1) / 2
        )]
        svPointColors <- svPointColors()

        # initialize the plot
        par(mar = c(4.1, 4.1, 0.1, 0.1))
        xlim <- c(min(targets$paddedStartI), max(targets$paddedEndI))
        ylim <- c(0, if(any(c("tt", "ta", "aa") %in% targetClasses())){
            targets[, max(paddedEndI)]
        } else {
            targets[, .(x = endI - paddedStartI + 1), by = regionName][, max(x)]
        })
        locationsPlot$initializeFrame(
            xlim = xlim,
            ylim = ylim,
            xlab = "SV Center (Mbp)",
            ylab = "SV Size (bp)",
            xaxs = "i", 
            yaxs = "i",
            xaxt = "n"              
        )
        # shade and demarcate the capture target regions
        targets[, {
            xinc <- ylim[2] / 2
            polygon(
                x = c(startI, startI + xinc, endI + xinc, endI, startI), 
                y = c(ylim[1], ylim[2], ylim[2], ylim[1], ylim[1]), 
                border = NA, col = "grey90"
            )
            polygon(
                x = c(startI, startI - xinc, endI - xinc, endI, startI), 
                y = c(ylim[1], ylim[2], ylim[2], ylim[1], ylim[1]), 
                border = NA, col = "grey90"
            )
            mtext(paste(regionName, chrom, sep = ","), side = 1, line = 2, at = centerI, cex = 1)
            paddedSize <- paddedEndI - paddedStartI + 1
            unit <- paddedSize / 10
            at <- seq(unit, paddedSize - unit, unit)
            axis(1, at = paddedStartI + at, labels = round((start - (startI - paddedStartI) + at) / 1e6, 2))
        }, by = regionName]    
        for(i in unlist(targets[, .(paddedStartI, startI, endI, paddedEndI)])){
            abline(-i * 2,  2, col = "grey60")
            abline( i * 2, -2, col = "grey60")
        }

        # plot the SV points on top
        locationsPlot$addPoints(
            x = svs[, centerI], 
            y = svs[, size], 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svPointColors$colors
        )

        # finally, plot any external comparator SVs on top of everything
        externalSvs <- externalSvs()
        if(!is.null(externalSvs)){     
            externalSvs <- externalSvs[chrom %in% targets$chrom]
            targets[, center := centerI + (start - startI)]
            targetIs <- mapply(function(svChrom, svCenter){ # find the target region closest to this CNV
                targets[, which.min(ifelse(chrom == svChrom, abs(center - svCenter), 1e9))]
            }, externalSvs$chrom, externalSvs$center)
            externalSvs[, centerI := center - targets[targetIs[.I], start - startI]] # convert to I coordinates
            locationsPlot$addPoints(
                x = externalSvs[, centerI], 
                y = externalSvs[, size], 
                pch = 1, 
                cex = 2, # stepSettings$Point_Size$value,
                col = CONSTANTS$plotlyColors$black
            ) 
        }
        selectedSv <- selectedSv()
        if(!is.null(selectedSv)){
            selectedSv[, ':='(
                size = abs(TARGET_POS_2 - TARGET_POS_1 + 1),
                centerI = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1 + 1) / 2
            )]
            locationsPlot$addPoints(
                x = selectedSv[, centerI], 
                y = selectedSv[, size], 
                pch = 1, 
                cex = 2, # stepSettings$Point_Size$value,
                col = CONSTANTS$plotlyColors$red
            )
        }

        # add a legend
        pointColorLegend(stepSettings, locationsPlot$settings, svPointColors)
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
nodesPlot <- junctionNodesPlotServer(svMols)
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
        nodesPlot$settings$replace(bm$outcomes$nodesPlotSettings)
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
        nodesPlotSettings = nodesPlot$settings$all_(),
        alignmentSettings = alignmentSettings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
