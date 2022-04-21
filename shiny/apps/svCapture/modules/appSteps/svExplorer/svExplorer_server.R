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
    template = c(file.path(app$sources$suiteGlobalDir, "sv_filters.yml"), id),
    fade = FALSE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues() # logical failure vectors keyed as [[sampleSet]]

#----------------------------------------------------------------------
# generate the list of all filtered SVs from all selected samples
#----------------------------------------------------------------------
targetClasses <- reactive({ SVX$targetClasses[[settings$SV_Filters()$Target_Class$value]] })
workingSvs <- reactive({
    getWorkingSvs(settings, sampleSelector)
})

#----------------------------------------------------------------------
# set SV plot point colors
#----------------------------------------------------------------------
svPointColors <- reactive({
    svs <- workingSvs()
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
pointColorLegend <- function(stepSettings, plotSettings, svPointColors){
    if(is.null(svPointColors$label)) return()
    legend(
        plotSettings$get('Plot_Frame', 'Legend_Placement'),
        svPointColors$label,
        pch = 20, 
        pt.cex = stepSettings$Point_Size$value * 1.25,
        col = svPointColors$color
    )
}

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
#     Classes 'data.table' and 'data.frame':  2 obs. of  11 variables:
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
        svs <- workingSvs()[, .(TARGET_POS_1, TARGET_POS_2)]
        svs[, ':='(
            size = abs(TARGET_POS_2 - TARGET_POS_1 + 1),
            center = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1 + 1) / 2
        )]
        svPointColors <- svPointColors()

        # initialize the plot
        par(mar = c(4.1, 4.1, 0.1, 0.1))
        plot(
            NA, NA, typ = "n",
            xlim = c(min(targets$paddedStartI), max(targets$paddedEndI)),
            ylim = c(0, if(any(c("tt", "ta", "aa") %in% targetClasses())){
                targets[, max(paddedEndI)]
            } else {
                targets[, .(x = endI - paddedStartI + 1), by = regionName][, max(x)]
            }),
            xlab = "SV Center (Mbp)",
            ylab = "SV Size (bp)",
            xaxs = "i", 
            yaxs = "i",
            xaxt = "n"
        )

        # shade and demarcate the capture target regions
        targets[, {
            polygon(
                x = c(startI, centerI, endI, startI), 
                y = c(0, size, 0, 0), 
                border = NA, col = "grey50"
            )
            mtext(paste(regionName, chrom, sep = ","), side = 1, line = 2.25, at = centerI, cex = 1)
            paddedSize <- paddedEndI - paddedStartI + 1
            unit <- paddedSize / 10
            at <- seq(unit, paddedSize - unit, unit)
            axis(1, at = paddedStartI + at, labels = round((start - (startI - paddedStartI) + at) / 1e6, 2))
        }, by = regionName]    
        for(i in unlist(targets[, .(paddedStartI, startI, endI, paddedEndI)])){
            abline(-i * 2,  2, col = "grey50")
            abline( i * 2, -2, col = "grey50")
        }

        # plot the SV points on top
        points(
            svs[, center], 
            svs[, size], 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svPointColors$colors
        )
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
        svs <- workingSvs()[, .(
            plotted = JXN_BASES != "*" & # MICROHOM_LEN meaningless if not a sequenced junction
                    JXN_TYPE  != "T",  # SV_SIZE meaningless for translocations
            MICROHOM_LEN, 
            SV_SIZE
        )]
        svPointColors <- svPointColors()
        par(mar = c(4.1, 4.1, 0.1, 1.1))
        plot(
            NA, 
            NA,
            typ = "n",
            xlim = c(-40, 40),
            ylim = log10(c(max(filters$Min_SV_Size$value, 1), filters$Max_SV_Size$value)),
            xlab = "Microhomology Length (bp)",
            ylab = "log10 SV Size (bp)"
        )
        abline(v = 0)
        points(
            jitter(svs[plotted == TRUE, MICROHOM_LEN]), 
            log10( svs[plotted == TRUE, SV_SIZE]), 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svPointColors$colors[svs[, plotted]]
        )
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
        workingSvs()[, .(
            svId = SV_ID,
            #---------------
            type = SVX$jxnTypes[JXN_TYPE, name],
            class = TARGET_CLASS,
            target = TARGET_REGION,
            #---------------
            nSmp = N_SAMPLES,
            nTot = N_TOTAL,
            nJxn = N_SPLITS,            
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
            uHom = MICROHOM_LEN,
            jxnSeq = JXN_BASES,
            svSize = SV_SIZE
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
# define bookmarking actions
# ----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    sampleSelector$setSampleSet(bm$input[['sampleSelector-sampleSet']]) 
    if(!is.null(bm$outcomes)) {
        outcomes <<- listToReactiveValues(bm$outcomes)
        sampleSelector$setSelectedSamples(bm$outcomes$samples)
        locationsPlot$settings$replace(bm$outcomes$locationsPlotSettings)
        propertiesPlot$settings$replace(bm$outcomes$propertiesPlotSettings)
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
        propertiesPlotSettings = propertiesPlot$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
