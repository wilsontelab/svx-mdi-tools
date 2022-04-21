#----------------------------------------------------------------------
# reactive components to show plots with QC values and distributions
# for all samples in a sampleSet, in order to explore data quality
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
libraryQCServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'libraryQC' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id
)
sampleSet <- sampleSetServer( # selectors to pick a sample set
    id = 'sampleSet',
    parentId = id
)
outcomes <- reactiveValues() # logical failure vectors keyed as [[sampleSet]]

#----------------------------------------------------------------------
# cascade from selected dataSource to its QC files
#----------------------------------------------------------------------
sourceIds <- reactive({
    req(sampleSet$assignments())
    unique(sampleSet$assignments()$Source_ID)
})
libraryMetrics <- reactive({
    req(sourceIds())
    startSpinner(session, 'libraryMetrics')   
    dt <- do.call(rbind, lapply(sourceIds(), function(sourceId){
        manifestFile <- getSourceFilePath(sourceId, "manifestFile")
        x <- fread(manifestFile)
        x[, ':='(
            Source_ID = sourceId,
            Sample_Unique_ID = paste(Project, Sample_ID, sep = ":"),
            efficiency = nSourceMolecules / nReadPairs
        )]
        x[, Sample_Name := getSampleNames(sampleUniqueIds = Sample_Unique_ID)]
        x
    }))
    stopSpinner(session, 'libraryMetrics')
    dt
})

#----------------------------------------------------------------------
# cascade to determine which sample are marked as failed, with initial defaults
#----------------------------------------------------------------------
failedStatus <- reactive({
    libraryMetrics <- libraryMetrics()
    req(libraryMetrics) 
    ss <- sampleSet$input$sampleSet
    if(is.null(outcomes[[ss]])){
        outcomes[[ss]] <- rep(FALSE, nrow(libraryMetrics))
        names(outcomes[[ss]]) <- libraryMetrics$Sample_Unique_ID
    }
    libraryMetrics()[, failed := outcomes[[ss]]]
    outcomes[[ss]]
})
output$nFailedLibraries <- renderText({
    failed <- unlist(failedStatus())
    req(!is.null(failed) && length(failed) > 0)
    paste(sum(failed),    'of', 
          length(failed), 'libraries marked as failed')
})

# #----------------------------------------------------------------------
# # systematically set and clear failure outcomes over many samples
# #----------------------------------------------------------------------
# observeEvent(input$applyFailedFilters, { # enforce filter settings _on top of_ existing failure marks
#     libraryMetrics <- libraryMetrics()
#     req(libraryMetrics) 
#     sourceId <- sampleSet$input$dataSource
#     failed <- failedStatus()
#     qct <- settings$QC_Thresholds()
#     failedReadCount <- libraryMetrics[, grouped   < qct$Minimum_Unique_Reads_M$value * 1e6]
#     failedAlignRate <- libraryMetrics[, alignRate < qct$Minimum_Alignment_Rate$value]
#     failedDupRate   <- libraryMetrics[, dupRate   > qct$Maximum_Duplication_Rate$value]
#     outcomes[[sourceId]] <- failed | failedReadCount | failedAlignRate | failedDupRate
#     names(outcomes[[sourceId]]) <- names(failed)
#     invalidateTable( invalidateTable() + 1 )
# })
# observeEvent(input$clearFailedMarks, { # clear _all_ failure marks, whether from settings or user click
#     libraryMetrics <- libraryMetrics()
#     req(libraryMetrics) 
#     sourceId <- sampleSet$input$dataSource
#     outcomes[[sourceId]] <- rep(FALSE, nrow(libraryMetrics))
#     invalidateTable( invalidateTable() + 1 )
# })

#----------------------------------------------------------------------
# stacked barplots of various library QC metrics
#----------------------------------------------------------------------
libraryMetricTypes <- list(
    nReadPairs  = "# Input Read Pairs",
    nSourceMolecules  = "# Source Molecules",
    onTargetCoverage  = "On-Target Coverage",
    offTargetCoverage  = "Off-Target Coverage",
    enrichment  = "Capture Enrichment",
    efficiency  = "Library Efficiency"
    # ,
    # nSvs  = "# SV Junctions"
)
libraryMetricsPlotData <- reactive({
    startSpinner(session, 'libraryMetricsPlotData')
    req(length(failedStatus()) > 0) 
    libraryMetrics <- libraryMetrics()[failed == FALSE]
    d <- lapply(names(libraryMetricTypes), function(col){
        d <- libraryMetrics[[col]]
        names(d) <- libraryMetrics$Sample_Name
        d
    })
    names(d) <- libraryMetricTypes
    stopSpinner(session, 'libraryMetricsPlotData')
    d
})  
getPlotRange <- function(d, i){
    max <- min(max(d), median(d) * settings$Plot_Limits()$Max_Fold_Median$value)
    c(0, max)
}
interactiveBarplotServer(
    'libraryMetricsPlot',
    libraryMetricsPlotData,
    orientation = 'vertical',
    shareAxis = list(x = TRUE),
    shareMargin = 0.01,
    lines = 'median',
    lineWidth = 1.5,
    range = getPlotRange
)

#----------------------------------------------------------------------
# scatter plot to explore metric relationships
#----------------------------------------------------------------------
metricRelationshipsData <- reactive({
    req(libraryMetrics()) 
    req(length(failedStatus()) > 0) 
    d <- libraryMetrics()[, .SD, .SDcols = c(input$xAxisMetric, input$yAxisMetric, 'failed')]
    names(d) <- c('x', 'y', 'failed')
    d
})
metricRelationshipsOverplotData <- reactive({ # nolint
    req(metricRelationshipsData())
    metricRelationshipsData()[failed == TRUE]
})
interactiveScatterplotServer(
    'metricRelationshipsPlot',
    metricRelationshipsData, 
    xtitle = reactive({ input$xAxisMetric }),
    ytitle = reactive({ input$yAxisMetric }),
    pointSize = 5, 
    overplot = metricRelationshipsOverplotData,
    overplotColor = "red",
    overplotPointSize = 5
)

# ----------------------------------------------------------------------
# summary table of all samples with FAILURE flag checkboxes
# ----------------------------------------------------------------------
invalidateTable <- reactiveVal(0)
librariesTable <- bufferedTableServer(
    id = 'librariesTable',
    parentId = id,
    parentInput = input,
    selection = 'single',
    tableData = reactive({ # construct the table data
        invalidateTable() # respond to table-external failure actions    
        startSpinner(session, 'bufferedTableServer')
        SDcols <- c('Sample_Name', names(libraryMetricTypes))
        lm <- libraryMetrics()[, .SD, .SDcols = c(SDcols, 'failed')]
        lm[, ':='(
            FAILED = tableCheckboxes(ns('libraryFailed'), failed ),
            Failed = failed,
            nReadPairs = commify(nReadPairs),
            nSourceMolecules = commify(nSourceMolecules),
            onTargetCoverage = commify(round(onTargetCoverage, 0)),
            offTargetCoverage = round(offTargetCoverage, 3),
            enrichment = commify(round(enrichment, 0)),
            efficiency = round(efficiency, 3)
            # ,
            # nSvs = commify(nSvs)
        )]
        stopSpinner(session, 'bufferedTableServer')
        lm[, .SD, .SDcols = c('FAILED', 'Failed', SDcols)]
    }),
    editBoxes = list( # handle the libraryFailed checkbox
        libraryFailed = list(
            type = 'checkbox',
            boxColumn = 1,
            rawColumn = 2,
            handler = function(checked){ # enter the new failure value into our outcomes
                ss <- sampleSet$input$sampleSet
                failed <- outcomes[[ss]]
                failed[checked$selectedRow] <- checked$newValue
                outcomes[[ss]] <- failed # must replace the entire value for outcomes invalidation
                checked
            }
        )
    )
)

#----------------------------------------------------------------------
# scatter traces of library insert size and family size distributions
#----------------------------------------------------------------------
getDistributionFilePath <- function(row, type){
    libraryMetrics <- libraryMetrics()
    req(libraryMetrics)
    fileName <- paste(libraryMetrics[row, Sample_ID], type, 'png', sep = ".")
    filePath <- expandSourceFilePath(libraryMetrics[row, Source_ID], fileName)
    if(!file.exists(filePath)){
        distributionsZip <- getSourceFilePath(libraryMetrics[row, Source_ID], 'distributionsZip')
        unzip(distributionsZip, exdir = dirname(distributionsZip))
    }
    filePath 
}
output$insertSizesPlot <- renderImage(
    {
        row <- librariesTable$rows_selected()
        req(row)
        list(
            src = getDistributionFilePath(row, 'insertSizes'),
            width = "100%"
        )
    }, 
    deleteFile = FALSE
)
output$familySizesPlot <- renderImage(
    {
        row <- librariesTable$rows_selected()
        req(row)
        list(
            src = getDistributionFilePath(row, 'family_sizes'),
            width = "100%"
        )
    }, 
    deleteFile = FALSE
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
