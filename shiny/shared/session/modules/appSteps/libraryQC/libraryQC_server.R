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
sampleSet <- sampleSetSourceServer( # selectors to pick a sample set
    id = 'sampleSet',
    parentId = id
)
outcomes <- reactiveValues() # outcomes[[sourceId]][i] <- TRUE if library i failed

# #----------------------------------------------------------------------
# # cascade from selected dataSource to its QC files
# #----------------------------------------------------------------------
# libraryMetrics <- reactive({
#     startSpinner(session, 'libraryMetrics')
#     sourceId <- sampleSet$input$dataSource
#     req(sourceId)
#     file <- getSourceFilePath(sourceId, 'sampleSummary')
#     dt <- fread(file)
#     sn1 <- getSampleNames(sampleIds = dt$sampleName) # hacky way to deal with inconsistent sample/library names in manifests # nolint
#     sn2 <- getSampleNames(sampleIds = dt$libraryName)
#     dt$insertSizesKey <- dt$sampleName
#     if(sum(is.na(sn1)) > sum(is.na(sn2))){
#         dt$Sample_ID <- dt$libraryName
#         dt$sampleName <- sn2
#     } else {
#         dt$Sample_ID <- dt$sampleName
#         dt$sampleName <- sn1
#     }
#     dt$efficiency <- dt$grouped / dt$readPairs
#     stopSpinner(session, 'libraryMetrics')
#     dt
# })
# insertSizes <- reactive({
#     startSpinner(session, 'insertSizes')
#     sourceId <- sampleSet$input$dataSource
#     req(sourceId)
#     file <- getSourceFilePath(sourceId, 'insertSizes')
#     req(file)
#     x <- fread(file)
#     stopSpinner(session, 'insertSizes')
#     x
# })

# #----------------------------------------------------------------------
# # cascade to determine which sample are marked as failed, with initial defaults
# #----------------------------------------------------------------------
# failedStatus <- reactive({
#     libraryMetrics <- libraryMetrics()
#     req(libraryMetrics)  
#     sourceId <- sampleSet$input$dataSource
#     if(is.null(outcomes[[sourceId]])) outcomes[[sourceId]] <- rep(FALSE, nrow(libraryMetrics))
#     names(outcomes[[sourceId]]) <- libraryMetrics$Sample_ID
#     outcomes[[sourceId]]
# })
# output$nFailedLibraries <- renderText({
#     failed <- failedStatus()
#     req(!is.null(failed) && length(failed) > 0)
#     paste(sum(failed),    'of', 
#           length(failed), 'libraries marked as failed')
# })

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

# #----------------------------------------------------------------------
# # stacked barplots of various library QC metrics
# #----------------------------------------------------------------------
# libraryMetricTypes <- list(
#     readPairs  = "# Input Read Pairs",
#     filtered   = "# Filtered Read Pairs",
#     grouped    = "# Unique Read Pairs",    
#     alignRate  = "Alignment Rate",
#     dupRate    = "Duplication Rate",
#     efficiency = "Library Efficiency"
# )
# isReadPairCount <- endsWith(unlist(libraryMetricTypes), "Read Pairs")
# libraryMetricsPlotData <- reactive({
#     startSpinner(session, 'libraryMetricsPlotData')
#     failed <- failedStatus() 
#     req(length(failed) > 0) 
#     libraryMetrics <- libraryMetrics()
#     d <- lapply(names(libraryMetricTypes), function(col){
#         d <- libraryMetrics[[col]]
#         names(d) <- libraryMetrics$libraryName
#         d[!failed]  
#     })
#     names(d) <- libraryMetricTypes
#     stopSpinner(session, 'libraryMetricsPlotData')
#     d
# })  
# getPlotRange <- function(d, i){
#     if(!isReadPairCount[i]) return( c(0, 1) ) # fractional quality metrics
#     max <- min(max(d), median(d) * settings$Plot_Limits()$Read_Count_Fold_Median$value)
#     c(0, max)
# }
# interactiveBarplotServer(
#     'libraryMetricsPlot',
#     libraryMetricsPlotData,
#     orientation = 'vertical',
#     shareAxis = list(x = TRUE),
#     shareMargin = 0.01,
#     lines = 'median',
#     lineWidth = 1.5,
#     range = getPlotRange
# )

# #----------------------------------------------------------------------
# # scatter traces of library insert sizes
# #----------------------------------------------------------------------
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

# #----------------------------------------------------------------------
# # summary table of all samples with FAILURE flag checkboxes
# #----------------------------------------------------------------------
# invalidateTable <- reactiveVal(0)
# bufferedTableServer(
#     id = 'librariesTable',
#     parentId = id,
#     parentInput = input,
#     selection = 'none',
#     tableData = reactive({ # construct the table data
#         startSpinner(session, 'bufferedTableServer')
#         SDcols <- c('libraryName', 'sampleName', names(libraryMetricTypes))
#         lm <- libraryMetrics()[, .SD, .SDcols = SDcols]
#         invalidateTable() # respond to table-external failure actions
#         lm[, ':='(
#             FAILED = tableCheckboxes(ns('libraryFailed'), isolate({ failedStatus() }) ),
#             Failed = isolate({ failedStatus() }),
#             alignRate  = round(alignRate, 3),
#             dupRate    = round(dupRate, 3),
#             efficiency = round(efficiency, 3)
#         )]
#         stopSpinner(session, 'bufferedTableServer')
#         lm[, .SD, .SDcols = c('FAILED', 'Failed', SDcols)]
#     }),
#     editBoxes = list( # handle the libraryFailed checkbox
#         libraryFailed = list(
#             type = 'checkbox',
#             boxColumn = 1,
#             rawColumn = 2,
#             handler = function(checked){ # enter the new failure value into our outcomes
#                 sourceId <- sampleSet$input$dataSource
#                 failed <- outcomes[[sourceId]]
#                 failed[checked$selectedRow] <- checked$newValue
#                 outcomes[[sourceId]] <- failed # must replace the entire value for outcomes invalidation
#                 checked
#             }
#         )
#     )
# )

# #----------------------------------------------------------------------
# # scatter plot to explore metric relationships
# #----------------------------------------------------------------------
# metricRelationshipsData <- reactive({
#     req(libraryMetrics()) 
#     d <- libraryMetrics()[, .SD, .SDcols = c(input$xAxisMetric, input$yAxisMetric)]
#     names(d) <- c('x', 'y')
#     d$failed <- failedStatus()
#     d
# })
# metricRelationshipsOverplotData <- reactive({ # nolint
#     req(metricRelationshipsData())
#     metricRelationshipsData()[failed == TRUE]
# })
# interactiveScatterplotServer(
#     'metricRelationshipsPlot',
#     metricRelationshipsData, 
#     xtitle = reactive({ input$xAxisMetric }),
#     ytitle = reactive({ input$yAxisMetric }),
#     pointSize = 5, 
#     overplot = metricRelationshipsOverplotData,
#     overplotColor = "red",
#     overplotPointSize = 5
# )

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
# observe({
#     bm <- getModuleBookmark(id, module, bookmark, locks)
#     req(bm)
#     settings$replace(bm$settings)
#     if(!is.null(bm$outcomes)) outcomes <<- listToReactiveValues(bm$outcomes)
# })

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

 #$ sampleN    : int  1 2 4 5 6 7 8 9 10 11 ...
 #$ sampleName : chr  "Sample_2251-LK-1" "Sample_2251-LK-2" "Sample_2251-LK-4" "Sample_2251-LK-5" ...
 #$ libraryName: chr  "2251-LK-1" "2251-LK-2" "2251-LK-4" "2251-LK-5" ...
 #$ cellName   : chr  "A1" "A2" "A4" "A5" ...
 #$ readPairs  : int  42380054 39037083 43699421 52921618 39102952 57393068 35436119 26553676 38667390 27447735 ...
 #$ filtered   : int  25704872 37102069 41310818 25952192 36460833 51499848 32529602 15185875 36615012 18491883 ...
 #$ grouped    : int  9236255 30787214 15085938 3837137 29360292 39219556 25104300 1090324 28183620 2303900 ...
 #$ alignRate  : num  0.607 0.95 0.945 0.49 0.932 ...
 #$ dupRate    : num  0.641 0.17 0.635 0.852 0.195 ...
 
 #$ sampleN    : int  1 2 3 4 5 6 7 8 9 10 ...
 #$ sampleName : chr  "Sample_1443-LK-1" "Sample_1443-LK-2" "Sample_1443-LK-3" "Sample_1443-LK-4" ...
 #$ libraryName: chr  "1443-LK-1" "1443-LK-2" "1443-LK-3" "1443-LK-4" ...
 #$ cellName   : logi  NA NA NA NA NA NA ...
 #$ readPairs  : int  7791719 9213974 7101140 8437780 8552801 9087751 7574990 8803321 12977043 10904835 ...
 #$ filtered   : int  344641 7108399 6407760 6892844 6818714 7397255 5895197 7049865 10341647 7892697 ...
 #$ grouped    : int  187900 6831362 6177410 6664854 6595123 7144605 5688860 6793934 9920731 7578402 ...
 #$ alignRate  : num  0.0442 0.7715 0.9024 0.8169 0.7972 ...
 #$ dupRate    : num  0.4548 0.039 0.0359 0.0331 0.0328 ...

 #$ insSize          : int  1 2 3 4 5 6 7 8 9 10 ...
 #$ Sample_1443-LK-10: num  0.00 5.28e-07 1.72e-06 1.85e-06 1.85e-06 ...
 #$ Sample_1443-LK-11: num  0.00 2.56e-07 2.56e-07 6.41e-07 2.56e-07 ...
