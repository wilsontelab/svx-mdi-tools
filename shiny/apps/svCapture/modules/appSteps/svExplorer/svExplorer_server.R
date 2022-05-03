#----------------------------------------------------------------------
# reactive components to filter and examine SV locations and junction sequences
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svExplorerServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'svExplorer' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    templates = list(
        file.path(app$sources$suiteGlobalDir, "sv_filters.yml"), 
        id,
        file.path(app$sources$suiteGlobalDir, "external_svs.yml")
    ),
    fade = FALSE
)
targetClasses <- reactive({ 
    SVX$targetClasses[[settings$SV_Filters()$Target_Class$value]] 
})
moduleOptions <- stepModuleInfo[[ app$config$appSteps[[id]]$module ]]
mapSettings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'mapSettings',
    parentId = id,
    templates = list(moduleOptions$mapSettings),
    fade = FALSE,
    title = "Base Map Settings",
    immediate = TRUE
)
alignmentSettings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'alignmentSettings',
    parentId = id,
    templates = list(moduleOptions$alignmentSettings),
    fade = FALSE,
    title = "Junction Alignment Settings",
    immediate = TRUE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
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

#----------------------------------------------------------------------
# load a set of external comparator SVs
#----------------------------------------------------------------------
externalSvs <- reactive({
    bedFile <- settings$External_SVs()$SVs_Bed_File$value
    req(bedFile)
    bedFile <- file.path(serverEnv$UPLOADS_DIR, bedFile)
    req(file.exists(bedFile))
    x <- fread(bedFile)[, 1:3]
    names(x) <- c("chrom", "start", "end")
    x[, ":="(
        size = end - start,
        center = start + (end - start) / 2
    )]
    x
})

#----------------------------------------------------------------------
# set SV plot point colors
#----------------------------------------------------------------------
svPointColors <- reactive({
    svs <- filteredSvs()
    req(svs)
    filters <- settings$SV_Filters()
    ps <- settings$Plot_Settings()
    setkey(SVX$jxnTypes, code)
    switch(ps$Point_Color$value,
        type = list(
            colors = svs[, SVX$jxnTypes[JXN_TYPE, color]],
            color = SVX$jxnTypes[!is.na(color) & name %in% filters$SV_Type$value, color],
            label = SVX$jxnTypes[!is.na(color) & name %in% filters$SV_Type$value, name]
        ),
        sample = {
            samples <- svs[N_SAMPLES == 1, sort(unique(SAMPLES))]
            sampleColors <- seq_along(samples)
            names(sampleColors) <- samples
            list(
                colors = svs[, ifelse(N_SAMPLES == 1, sampleColors[SAMPLES], "grey")],
                color = sampleColors,
                label = samples
            )
        },
        target = {
            targets <- svs[!grepl(',', TARGET_REGION), sort(unique(TARGET_REGION))]
            targetColors <- seq_along(targets)
            names(targetColors) <- targets
            list(
                colors = svs[, ifelse(grepl(',', TARGET_REGION), "grey", targetColors[TARGET_REGION])],
                color = targetColors,
                label = targets
            )
        }, 
        duplex = list(
            colors = svs[, ifelse(N_DUPLEX_GS > 0, "green3",  "red3")],
            color = c("green3",  "red3"),
            label = c("Duplex", "Single Strand")
        ),
        blue = list(
            colors = "blue", 
            color = NULL,
            label = NULL
        )
    )
})
pointColorLegend <- function(stepSettings, plotSettings, svPointColors, exclude = character()){
    if(is.null(svPointColors$label)) return()
    is <- which(!(svPointColors$label %in% exclude))
    legend(
        plotSettings$get('Plot_Frame', 'Legend_Placement'),
        svPointColors$label[is],
        pch = 20, 
        pt.cex = stepSettings$Point_Size$value * 1.25,
        col = svPointColors$color[is]
    )
}

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
locationsPlot <- staticPlotBoxServer(
    'svLocations', 
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
            center = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1 + 1) / 2
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
        plot(
            NA, NA, typ = "n",
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
            mtext(paste(regionName, chrom, sep = ","), side = 1, line = 2.25, at = centerI, cex = 1)
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
        points(
            svs[, center], 
            svs[, size], 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svPointColors$colors
        )

        # finally, plot any external comparator SVs on top of everything
        externalSvs <- externalSvs()
        if(!is.null(externalSvs)){ 
# Classes 'data.table' and 'data.frame':  16 obs. of  5 variables:
#  $ chrom : chr  "chr10" "chr10" "chr1" "chr1" ...
#  $ start : int  51616453 51628463 71643114 71830439 51508820 78602605 71861295 78488898 71718794 51566939 ...
#  $ end   : int  51648536 51647510 71768387 72047181 51744414 78946165 71886554 78576240 72164112 51695648 ...
#  $ size  : int  32083 19047 125273 216742 235594 343560 25259 87342 445318 128709 ...
#  $ center: num  51632495 51637987 71705751 71938810 51626617 ...
#  - attr(*, ".internal.selfref")=<externalptr>
# Classes 'data.table' and 'data.frame':  2 obs. of  11 variables:
#  $ chrom       : chr  "chr3" "chr16"
#  $ start       : int  60300384 78456040
#  $ end         : int  60700384 78856040
#  $ regionName  : chr  "FHIT" "WWOX"
#  $ chromI      : int  3 16
#  $ paddedStartI: int  1 2000001
#  $ startI      : int  800001 2800001
#  $ centerI     : int  1000000 3000000
#  $ endI        : int  1200000 3200000
#  $ paddedEndI  : int  2000000 4000000
#  $ size        : int  400000 400000          
            externalSvs <- externalSvs[chrom %in% targets$chrom]
            # externalSvs <- externalSvs[, {
            #     svChrom <- chrom
            #     targets[targets$chrom == ]
            # }, by = chrom]

            points(
                externalSvs[, center], 
                externalSvs[, size], 
                pch = "X"
                # , 
                # cex = stepSettings$Point_Size$value,
                # col = svPointColors$colors
            )  
        }

        # add a legend
        pointColorLegend(stepSettings, locationsPlot$settings, svPointColors)
    }
)

#----------------------------------------------------------------------
# SV junction properties plot
#----------------------------------------------------------------------
propertiesPlot <- staticPlotBoxServer(
    'svProperties', 
    legend = TRUE,
    immediate = TRUE,
    create = function(...){
        filters <- settings$SV_Filters()
        stepSettings <- settings$Plot_Settings()   
        svPointColors <- svPointColors()        
        svs <- filteredSvs()[, .(
            plotted = JXN_BASES != "*", # MICROHOM_LEN meaningless if not a sequenced junction
            color = svPointColors$colors,
            MICROHOM_LEN, 
            SV_SIZE,
            JXN_TYPE
        )]
        svs[JXN_TYPE == "T", SV_SIZE := rnorm(.N, filters$Max_SV_Size$value, filters$Max_SV_Size$value / 10)]
        par(mar = c(4.1, 4.1, 0.1, 1.1))
        maxX <- stepSettings$Max_Microhomology$value
        xlim <- c(-maxX, maxX)
        ylim <- log10(c(max(filters$Min_SV_Size$value, 1), filters$Max_SV_Size$value * 1.5))
        plot(
            NA, 
            NA,
            typ = "n",
            xlim = xlim,
            ylim = ylim,
            xlab = "Microhomology Length (bp)",
            ylab = "log10 SV Size (bp)"
        )
        abline(h = seq(0, 10, 1), col = "grey60")
        abline(v = seq(-100, 100, 10), col = "grey60")
        abline(v = 0)
        points(
            jitter(svs[plotted == TRUE, MICROHOM_LEN], amount = 0.5), # spread points from (i-0.5):(i+0.5)
            log10( svs[plotted == TRUE, SV_SIZE]), 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svs[plotted == TRUE, color]
        )

        # add text to denote microhomology vs. insertions
        y <- ylim[1] + (ylim[2] - ylim[1]) * 0.035
        # text(xlim[1], y, "insertions", pos = 4, offset = 0.05, cex = 0.95)
        # text(xlim[2], y, expression(paste(mu, "homology")), pos = 2, offset = 0.05, cex = 0.95)
        mtext("insertion", side = 1, line = 2, at = xlim[1], adj = 0, cex = 0.95)
        mtext(expression(paste(mu, "homology")), side = 1, line = 2, at = xlim[2], adj = 1, cex = 0.95)

        # add a legend
        pointColorLegend(stepSettings, propertiesPlot$settings, svPointColors)
    }
)

# ----------------------------------------------------------------------
# summary table of all filtered SVs from all selected samples
# ----------------------------------------------------------------------
invalidateTable <- reactiveVal(0)
svsTable <- bufferedTableServer(
    id = 'svsTable',
    parentId = id,
    parentInput = input,
    selection = 'single',
    tableData = reactive({
        invalidateTable()
        setkey(SVX$jxnTypes, code)
        filteredSvs()[, .(
            svId = SV_ID,
            #---------------
            type = SVX$jxnTypes[JXN_TYPE, name],
            class = TARGET_CLASS,
            target = TARGET_REGION,
            #---------------
            nSmp = N_SAMPLES,
            nTot = N_TOTAL,
            nSplit = N_SPLITS,            
            nGap = N_GAPS,
            nClip = N_OUTER_CLIPS,
            nDpx = N_DUPLEX_GS,
            #---------------
            count = STRAND_COUNT_GS,
            count1 = STRAND_COUNT1,
            count2 = STRAND_COUNT2,
            #---------------
            chr1 = CHROM_1,
            pos1 = POS_1,
            chr2 = CHROM_2,
            pos2 = POS_2,
            #---------------
            svSize = SV_SIZE,
            uHom = MICROHOM_LEN,
            jxnSeq = JXN_BASES
        )]
    })
    # 'SHARED_PROPER' =  'double',
    # 'SHARED_PROPER_GS'  =  'double',
    # #---------------
    # 'N_AMBIGUOUS'   =  'integer',
    # 'N_DOWNSAMPLED' =  'integer',
    # 'N_COLLAPSED'   =  'integer',
    #     lm[, ':='(
    #         FAILED = tableCheckboxes(ns('libraryFailed'), failed ),
    #         Failed = failed,
    #         nReadPairs = commify(nReadPairs),
    #         nSourceMolecules = commify(nSourceMolecules),
    #         onTargetCoverage = commify(round(onTargetCoverage, 0)),
    #         offTargetCoverage = round(offTargetCoverage, 3),
    #         enrichment = commify(round(enrichment, 0)),
    #         efficiency = round(efficiency, 3)
    #     )]
    # editBoxes = list( # handle the libraryFailed checkbox
    #     libraryFailed = list(
    #         type = 'checkbox',
    #         boxColumn = 1,
    #         rawColumn = 2,
    #         handler = function(checked){ # enter the new failure value into our outcomes
    #             ss <- sampleSet$input$sampleSet
    #             failed <- outcomes[[ss]]
    #             failed[checked$selectedRow] <- checked$newValue
    #             outcomes[[ss]] <- failed # must replace the entire value for outcomes invalidation
    #             checked
    #         }
    #     )
    # )
)

# ----------------------------------------------------------------------
# plot of outer endpoints of all molecules matching the selected junction
# ----------------------------------------------------------------------
junctionNodesPlot <- staticPlotBoxServer(
    'junctionNodes', 
    legend = TRUE,
    points = TRUE,
    margins = TRUE,
    immediate = TRUE,
    create = function(...){
        x <- svMols()
        req(x)
        xAllowedLim <- x$sv$POS_1 + c(-1, 1) * x$maxTLen
        yAllowedLim <- x$sv$POS_2 + c(-1, 1) * x$maxTLen
        x$mols <- x$mols[ # in case molecule has >1 jxn; only plot this junction
            OUT_POS1 %between% xAllowedLim & 
            OUT_POS2 %between% yAllowedLim
        ]
        junctionNodesPlot$initializeFrame(
            xlim = range(c(x$sv$POS_1, x$mols$OUT_POS1)) / 1e6,
            ylim = range(c(x$sv$POS_2, x$mols$OUT_POS2)) / 1e6,
            xlab = "Junction Coordinate 1 (Mbp)",
            ylab = "Junction Coordinate 2 (Mbp)"
        )
        abline(h = x$sv$POS_2 / 1e6) # crosshairs at the junction point
        abline(v = x$sv$POS_1 / 1e6)   
        junctionNodesPlot$addPoints(
            x = c(x$sv$POS_1, x$mols$OUT_POS1) / 1e6, 
            y = c(x$sv$POS_2, x$mols$OUT_POS2) / 1e6,
            col = c(
                CONSTANTS$plotlyColors$red, # red dot at the junction point
                SVX$getMolColors(x$mols$NODE_CLASS)
            )
        )
        junctionNodesPlot$addLegend(
            legend = names(SVX$nodeClassColors),
            col   = unlist(SVX$nodeClassColors)
        )
    }
)

# ----------------------------------------------------------------------
# colored-image representation of all molecules supporting and SV junction
# ----------------------------------------------------------------------
output$junctionMapImage <- renderImage({
    map <- junctionMap()
    req(map)
    pngFile <- file.path(sessionDirectory, "junctionMapImage.png")
    pixels <- mapSettings$get("Map_Settings", "Pixels_Per_Base")
    suppressWarnings(
        imager::as.cimg(map$image[map$usedPos, , ]) %>% 
        expandImg(h = pixels, v = pixels) %>%
        imager::save.image(pngFile)        
    )
    list(src = pngFile)
}, deleteFile = FALSE)

# ----------------------------------------------------------------------
# text representation of all an SV junction consensus sequence vs. genome references
# ----------------------------------------------------------------------
output$junctionAlignment <- renderText({
    getJunctionAlignment(
        junctionMap(),
        alignmentSettings$get("Alignment_Settings", "Bases_Per_Line"),
        alignmentSettings$get("Alignment_Settings", "Display_Mode")
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
