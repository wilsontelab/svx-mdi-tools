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
settings <- stepSettingsServer( # display settings not stored in the UI, exposed by gear icon click
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
    startSpinner(session, 'libraryMetrics')
# sampleSet$assignments
# Classes 'data.table' and 'data.frame':  6 obs. of  5 variables:
#  $ Source_ID: chr  "6741c8f40c01f6ed4760406898b1607c" "6741c8f40c01f6ed4760406898b1607c" "6741c8f40c01f6ed4760406898b1607c" "6741c8f40c01f6ed4760406898b1607c" ...
#  $ Project  : chr  "high_APH_040422" "high_APH_040422" "high_APH_040422" "high_APH_040422" ...
#  $ Sample_ID: chr  "HCT_0.2APH_2APH_M" "HCT_0.2APH_G2" "HCT_0.2APH__M" "HCT_2APH_M" ...
#  $ Category1: int  2 1 2 2 1 2
#  $ Category2: int  4 2 2 3 1 1
    unique(sampleSet$assignments()$Source_ID)
})
libraryMetrics <- reactive({
    req(sourceIds())
    dt <- do.call(rbind, lapply(sourceIds(), function(sourceId){
        manifestFile <- getSourceFilePath(sourceId, "manifestFile")
        x <- fread(manifestFile)
        x[, ':='(
            Source_ID = sourceId,
            Sample_Name = paste(Project, Sample_ID, sep = ":"),
            efficiency = nSourceMolecules / nReadPairs
        )]
        x
    }))
    stopSpinner(session, 'libraryMetrics')
    dt
})

observeEvent(libraryMetrics(), {
# libraryMetrics()
# Classes 'data.table' and 'data.frame':  6 obs. of  12 variables:
#  $ Project          : chr  "high_APH_040422" "high_APH_040422" "high_APH_040422" "high_APH_040422" ...
#  $ Sample_ID        : chr  "HCT_0.2APH_2APH_M" "HCT_0.2APH_G2" "HCT_0.2APH__M" "HCT_2APH_M" ...
#  $ Description      : chr  "HCT_0.2APH_2APH_M" "HCT_0.2APH_G2" "HCT_0.2APH__M" "HCT_2APH_M" ...
#  $ maxTLen          : int  790 844 716 762 840 843
#  $ nReadPairs       : int  134723261 103110289 173540492 162222559 182287552 167856108
#  $ nSourceMolecules : int  9018381 8339919 16452974 9553548 10312938 11943757
#  $ onTargetCoverage : num  2444 2275 2623 2540 2602 ...
#  $ offTargetCoverage: num  0.096 0.096 0.598 0.111 0.147 0.237
#  $ enrichment       : num  25462 23695 4386 22883 17702 ...
#  $ nSvs             : int  5053 3887 5822 4755 4378 4864
#  $ Yield            : int  134723261 103110289 173540492 162222559 182287552 167856108
#  $ sourceId         : chr  "6741c8f40c01f6ed4760406898b1607c" "6741c8f40c01f6ed4760406898b1607c" "6741c8f40c01f6ed4760406898b1607c" "6741c8f40c01f6ed4760406898b1607c" ...
    dstr(libraryMetrics())
})

# insertSizeFiles <- reactive({
#     req(sourceIds())
#     dt <- do.call(rbind, lapply(sourceIds(), function(sourceId){
#         distributionsZip <- getSourceFilePath(sourceId, "distributionsZip")
#         fread(manifestFile)
#     }))


#     dt <- do.call(rbind, lapply(sourceIds, function(sourceId){
#         manifestFile <- getSourceFilePath(sourceId, "manifestFile")
#         fread(manifestFile)
#     }))


#     sourceId <- sampleSet$input$dataSource
#     req(sourceId)
#     file <- getSourceFilePath(sourceId, 'insertSizes')
#     req(file)
#     x <- fread(file)
#     stopSpinner(session, 'insertSizes')
#     x
# })

#----------------------------------------------------------------------
# cascade to determine which sample are marked as failed, with initial defaults
#----------------------------------------------------------------------
failedStatus <- reactive({
    libraryMetrics <- libraryMetrics()
    req(libraryMetrics) 
    ss <- sampleSet$input$sampleSet
    if(is.null(outcomes[[ss]])){
        outcomes[[ss]] <- rep(FALSE, nrow(libraryMetrics))
        names(outcomes[[ss]]) <- libraryMetrics$Sample_Name
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
    efficiency  = "Library Efficiency",
    nSvs  = "# SV Junctions"
)

isReadPairCount <- endsWith(unlist(libraryMetricTypes), "Read Pairs")

libraryMetricsPlotData <- reactive({
    startSpinner(session, 'libraryMetricsPlotData')
    req(length(failedStatus()) > 0) 
    libraryMetrics <- libraryMetrics()[failed == FALSE]
    d <- lapply(names(libraryMetricTypes), function(col){
        d <- libraryMetrics[[col]]
        names(d) <- libraryMetrics$Sample_ID
        d
    })
    names(d) <- libraryMetricTypes
    stopSpinner(session, 'libraryMetricsPlotData')
    d
})  
getPlotRange <- function(d, i){
    # if(!is.integer(d)) return( c(0, 1) ) # fractional quality metrics
    # if(!isReadPairCount[i]) return( c(0, 1) ) # fractional quality metrics
    max <- min(max(d), median(d) * settings$Plot_Limits()$Read_Count_Fold_Median$value)
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
bufferedTableServer(
    id = 'librariesTable',
    parentId = id,
    parentInput = input,
    selection = 'none',
    tableData = reactive({ # construct the table data
        invalidateTable() # respond to table-external failure actions    
        startSpinner(session, 'bufferedTableServer')
        SDcols <- c('Sample_ID', names(libraryMetricTypes))
        lm <- libraryMetrics()[, .SD, .SDcols = c(SDcols, 'failed')]
        lm[, ':='(
            FAILED = tableCheckboxes(ns('libraryFailed'), failed ),
            Failed = failed,
            nReadPairs = commify(nReadPairs),
            nSourceMolecules = commify(nSourceMolecules),
            onTargetCoverage = commify(round(onTargetCoverage, 0)),
            offTargetCoverage = round(offTargetCoverage, 3),
            enrichment = commify(round(enrichment, 0)),
            efficiency = round(efficiency, 3),
            nSvs = commify(nSvs)
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



# insertSizesPlotData <- reactive({
#     startSpinner(session, 'insertSizesPlotData')
#     failed <- failedStatus() 
#     req(length(failed) > 0) 
#     libraryMetrics <- libraryMetrics()
#     insertSizes <- insertSizes()
#     req(insertSizes)
#     i <- 1
#     x <- do.call(rbind, lapply(libraryMetrics$insertSizesKey[!failed], function(sampleId){
#         dt <- data.table(
#             x = insertSizes$insSize,
#             y = insertSizes[[sampleId]]
#         )
#         if(i %% 2) dt <- dt[order(-x)]
#         i <<- i + 1
#         dt
#     })) 
#     stopSpinner(session, 'insertSizesPlotData')
#     x
# })
# insertSizeRange <- reactive({ # attempt to focus on the peak of the insert size distribution
#     d <- insertSizesPlotData()
#     d <- d[, .(y = sum(y)), keyby = x]
#     cum <- cumsum(d$y) / sum(d$y, na.rm = TRUE)
#     c(0, which.max(d$x[cum <= 0.99]))
# })
# interactiveScatterplotServer(
#     'insertSizesPlot',
#     insertSizesPlotData, 
#     xtitle = "InsertSize",
#     ytitle = "Frequency",
#     xrange = insertSizeRange,
#     yrange = range_pad,
#     mode = 'lines',
#     lineWidth = 1
# )

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
