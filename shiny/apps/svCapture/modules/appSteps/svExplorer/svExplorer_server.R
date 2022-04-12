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
    parentId = id
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
    tf  <- settings$Target_Filters()
    svf <- settings$SV_Filters()
    x <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
        svFile <- getSourceFilePath(sourceId, "structuralVariants")
        fillSvxCache(sourceId, 'svs', svFile) 

        # TODO: sample filters via selected assignments
        
        svxCache$svs[[svFile]][
            TARGET_CLASS %in% tf$Target_Class$value &
            JXN_TYPE %in% unlist(SVX$jxnTypeKey[svf$SV_Type$value]) &
            (JXN_TYPE == "T" | SV_SIZE >= svf$Min_SV_Size$value) &
            STRAND_COUNT >= svf$Min_Read_Count$value &
            N_SAMPLES >= svf$Min_Samples_With_SV$value &
            N_SAMPLES <= svf$Max_Samples_With_SV$value & 
            TRUE
            , .SD,
           .SDcols = names(SVX$find$structural_variants) 
        ]
    }))
    stopSpinner(session, 'workingSvs')
    x
})

output$svLocationsPlot <- renderPlot({
    svs <- workingSvs()[, .(TARGET_POS_1, TARGET_POS_2)]
    plot(svs$TARGET_POS_1, svs$TARGET_POS_2, pch = ".")
    abline(0, 1)
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
