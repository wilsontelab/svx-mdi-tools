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
settings <- stepSettingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
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
workingSvs <- reactive({
    assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
    req(assignments)
    req(nrow(assignments) > 0)
    startSpinner(session, 'workingSvs')
    filters <- settings$SV_Filters()
    samples <- sampleSelector$selectedSamples()
    isAllSamples <- all(sampleSelector$allSamples() %in% samples)
    x <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
        
        # load table into cache
        svFile <- fillSvxCache_svs(sourceId) 

        # apply SV filters
        x <- svxCache$svs[[svFile]][
            TARGET_CLASS %in% filters$Target_Class$value &
            JXN_TYPE %in% unlist(SVX$jxnTypeKey[filters$SV_Type$value]) &
            (JXN_TYPE == "T" | 
             (SV_SIZE >= filters$Min_SV_Size$value & 
              SV_SIZE <= filters$Max_SV_Size$value)) &
            STRAND_COUNT >= filters$Min_Read_Count$value &
            N_SAMPLES >= filters$Min_Samples_With_SV$value &
            N_SAMPLES <= filters$Max_Samples_With_SV$value
        ]

        # apply sample filters, unless showing all samples
        if(!isAllSamples){
            x[, matchesSamples := FALSE]
            project <- assignments[Source_ID == sourceId][1, Project]
            for(sample in samples){
                d <- strsplit(sample, ":")[[1]]
                if(d[1] == project) x[, 
                    matchesSamples := matchesSamples | x[[d[2]]] > 0
                ]
            }
        } else {
            x[, matchesSamples := TRUE]
        }
        x[matchesSamples == TRUE, .SD, .SDcols = names(SVX$find$structural_variants) ]
    }))
    stopSpinner(session, 'workingSvs')
    x
})

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
output$svLocationsPlot <- renderPlot({
    # lines(c(0, maxI), c(1, 1), col = "black", lwd = 2)
    # for(i in 1:nCaptureTargets){
    #     # diagonal lines
    #     i0 <- captureTargets[i, 'iOffset']
    #     i1 <- i0 + projectInfo$REGION_PADDING
    #     i2 <- i1 + captureTargets[i, 'length']
    #     i3 <- i2 + projectInfo$REGION_PADDING
    #     addGridline(i1, 1,  col = "grey60")
    #     addGridline(i1, -1, col = "grey60")
    #     addGridline(i2, 1,  col = "grey60")
    #     addGridline(i2, -1, col = "grey60")  
    #     addGridline(i0, 1)
    #     addGridline(i3, -1)    
        
    #     # coordinate axis per region
    #     unit <- 50000
    #     at <- seq(unit, captureTargets[i, 'paddedLength'] - unit, unit)
    #     axis(1, at = at + captureTargets[i, 'iOffset'], labels = round((captureTargets[i,'start'] + at) / 1e6, 2))
        
    #     # region name label
    #     mtext(paste(captureTargets[i, 'region'], captureTargets[i, 'chrom'], sep = ","),
    #           side = 1, line = 2.25, at = captureTargets[i, 'iOffset'] + captureTargets[i, 'paddedLength'] / 2,
    #           cex = regionLabelCex)
        
    #     # coordinate name
    #     mtext("SV center (Mbp)", side = 1, line = 3.5, at = maxI / 2, cex = axisLabelCex)
    # }

    # load targets bed, assumed to be the same for all sample sources
    assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
    req(assignments)
    req(nrow(assignments) > 0)
    sourceId <- assignments[, Source_ID[1]]
    targetsBed <- fillSvxCache_targets(sourceId)
    targets <- svxCache$targets[[targetsBed]]
    maxPaddedSize <- targets[, .(size = paddedEndI - paddedStartI), by = regionName][, max(size)]

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

    svs <- workingSvs()[, .(TARGET_POS_1, TARGET_POS_2)]
    svs[, ':='(
        size = abs(TARGET_POS_2 - TARGET_POS_1),
        center = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1) / 2
    )]
    plot(
        svs[, center], 
        svs[, size], 
        pch = 16, 
        cex = 0.25,
        xlim = c(min(targets$paddedStartI), max(targets$paddedEndI)),
        ylim = c(0, maxPaddedSize),
        xlab = "SV Center (bp)",
        ylab = "SV Size (bp)"
    )
    for(i in unlist(targets[, .(paddedStartI, startI, endI, paddedEndI)])){
        abline(-i * 2, 2)
        abline(i * 2, -2)
    }
})

#----------------------------------------------------------------------
# SV junction properties plot
#----------------------------------------------------------------------
output$svPropertiesPlot <- renderPlot({
    filters <- settings$SV_Filters()
    svs <- workingSvs()[
        JXN_BASES != "*" & # MICROHOM_LEN meaningless if not a sequenced junction
        JXN_TYPE != "T",   # SV_SIZE meaningless for translocations
        .(MICROHOM_LEN, SV_SIZE)
    ]
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
        jitter(svs$MICROHOM_LEN), 
        log10(svs$SV_SIZE), 
        pch = 16, 
        cex = 0.25,
        col = "blue"
    )
})

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
        workingSvs()[, .(
            svId = SV_ID,
            #---------------
            type = unlist(SVX$jxnTypeName[JXN_TYPE]),
            class = TARGET_CLASS,
            target = TARGET_REGION,
            #---------------
            nSmp = N_SAMPLES,
            nTot = N_TOTAL,
            nGap = N_GAPS,
            nJxn = N_SPLITS,
            nClip = paste0('(', N_OUTER_CLIPS, ')'),
            nDpx = paste0(N_DUPLEX_GS, '(', N_DUPLEX, ')'),
            count = paste0(STRAND_COUNT_GS, '(', STRAND_COUNT, ')'),


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
 
    # 'STRAND_COUNT1' =  'integer',
    # 'STRAND_COUNT2' =  'integer',
    # 'TARGET_CLASS'  =  'character',
    # 'SHARED_PROPER' =  'double',
    # 'SHARED_PROPER_GS'  =  'double',
    # #---------------
    # 'N_AMBIGUOUS'   =  'integer',
    # 'N_DOWNSAMPLED' =  'integer',
    # 'N_COLLAPSED'   =  'integer',
    # tableData = reactive({ # construct the table data
    #     invalidateTable() # respond to table-external failure actions    
    #     startSpinner(session, 'bufferedTableServer')
    #     SDcols <- c('Sample_Name', names(libraryMetricTypes))
    #     lm <- libraryMetrics()[, .SD, .SDcols = c(SDcols, 'failed')]
    #     lm[, ':='(
    #         FAILED = tableCheckboxes(ns('libraryFailed'), failed ),
    #         Failed = failed,
    #         nReadPairs = commify(nReadPairs),
    #         nSourceMolecules = commify(nSourceMolecules),
    #         onTargetCoverage = commify(round(onTargetCoverage, 0)),
    #         offTargetCoverage = round(offTargetCoverage, 3),
    #         enrichment = commify(round(enrichment, 0)),
    #         efficiency = round(efficiency, 3)
    #         # ,
    #         # nSvs = commify(nSvs)
    #     )]
    #     stopSpinner(session, 'bufferedTableServer')
    #     lm[, .SD, .SDcols = c('FAILED', 'Failed', SDcols)]
    # }),
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
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    samples  = sampleSelector$selectedSamples,
    # outcomes = reactive({ reactiveValuesToList(outcomes) }),
    outcomes = reactive({ list(
        samples = sampleSelector$selectedSamples()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
