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

# TODO: move cache to server level, not session
cache <- reactiveValues(svs = list(), molecules = list())

workingSvs <- reactive({
    assignments <- sampleSelector$selectedAssignments()
    req(assignments)
    req(nrow(assignments) > 0)
    startSpinner(session, 'workingSvs')
    activeSvFiles <- list()
    x <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
        svFile <- getSourceFilePath(sourceId, "structuralVariants")
        reportProgress(svFile)
        activeSvFiles[[svFile]] <<- TRUE        
        if(is.null(cache$svs[[svFile]])) cache$svs[[svFile]] <- fread(svFile)



        # TODO: implement SV filters via Settings
        # and sample filters via selected assignments

        cache$svs[[svFile]][uniqueId %in% ]
    }))
    
    for(svFile in names(cache$svs)){
        if(is.null(activeSvFiles[[svFile]])) cache$svs[[svFile]] <- NULL
    }
    stopSpinner(session, 'workingSvs')
    x
    # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
})

output$svLocationsPlot <- renderPlot({
    svs <- workingSvs()[JXN_TYPE == "T", .(TARGET_POS_1, TARGET_POS_2)]
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

    # TODO: parse the information required for the table
    tableData = workingSvs #sampleSelector$selectedAssignments

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
    if(!is.null(bm$outcomes)) outcomes <<- listToReactiveValues(bm$outcomes)
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = reactive({ reactiveValuesToList(outcomes) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
