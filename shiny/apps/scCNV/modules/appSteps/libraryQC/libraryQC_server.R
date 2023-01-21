#----------------------------------------------------------------------
# static components to show plots with QC values and distributions
# for all samples in a sampleSet, to explore run quality and reject cells
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
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
outcomes <- reactiveValues() # outcomes[[sourceId]][i] <- TRUE if library i failed

#----------------------------------------------------------------------
# cascade from selected dataSource to its QC data
#----------------------------------------------------------------------
colData <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading cell data")
    x <- loadSampleRds(sourceId)
    stopSpinner(session)
    x$colData
})

#----------------------------------------------------------------------
# make the QC plots
#----------------------------------------------------------------------
cellCol <- reactive({
    d <- colData()
    cols <- CONSTANTS$plotlyColors
    ifelse(d$bad,         cols$red, 
    ifelse(!d$keep,       cols$orange,
    ifelse(d$replicating, cols$blue,
                          cols$green)))
})
qcPlot <- function(plot, xcol, ycol, xlab, ylab, jitterX = FALSE, log10 = FALSE){
    d <- colData()
    x <- d[[xcol]]
    y <- d[[ycol]]
    if(log10) x <- log10(x)
    if(jitterX) x <- jitter(x)
    par(mar = c(4.1, 4.1, 0.1, 0.1))
    plot$initializeFrame(
        xlab = xlab,
        ylab = ylab,
        xlim = range(x, na.rm = TRUE),
        ylim = range(y, na.rm = TRUE)
    )
    plot$addPoints( # addLines follows the same pattern
        x = x,
        y = y,
        pch = 16,
        cex = 0.75,
        col = cellCol()
    )
}

# mdiInteractivePlotServer(
#     'nReads_vs_alignRate',     
#     # click = TRUE,
#     # brush = TRUE,
#     # delay = 500,
#     #------------------------------
#     contents = reactive({

#     }) # a reactive that returns:
#     # contents = reactive({ list(
#     #     pngFile = path, # OR plotArgs
#     #     plotArgs = list(
#     #          ... # first element must have a "plot.class" method and be named "x" or unnamed
#     #     ),
#     #     layout = list(
#     #         width = pixels,
#     #         height = pixels,
#     #         pointsize = integer, # defaults to 8
#     #         dpi = integer, # defaults to 96
#     #         mai = integer vector,
#     #         xlim = range, # OR can be read from plotArgs
#     #         ylim = range
#     #     ),
#     #     parseLayout = function(x, y) list(x, y, layout) # to convert to plot space in a multi-plot layout
#     # ) })
# )


arPlot <- staticPlotBoxServer(
    'nReads_vs_alignRate', 
    points = TRUE,
    create = function(){
        qcPlot(
            arPlot,
            xcol = "alignRate",
            ycol = "total_num_reads",
            xlab = "Alignment Rate",
            ylab = "Total # of Reads"
        )
        # abline(v = 0) # an example of direct plot manipulation
    }
)
drPlot <- staticPlotBoxServer(
    'nReads_vs_dupRate', 
    points = TRUE,
    create = function(){
        qcPlot(
            drPlot,
            xcol = "dupRate",
            ycol = "total_num_reads",
            xlab = "Duplication Rate",
            ylab = "Total # of Reads"
        )
        # abline(v = 0) # an example of direct plot manipulation
    }
)
wpPlot <- staticPlotBoxServer(
    'nReads_vs_windowPower', 
    points = TRUE,
    create = function(){
        qcPlot(
            arPlot,
            xcol = "windowPower",
            ycol = "total_num_reads",
            xlab = "Window Power",
            ylab = "Total # of Reads",
            jitterX = TRUE
        )
        # abline(v = 0) # an example of direct plot manipulation
    }
)
cnsdPlot <- staticPlotBoxServer(
    'nReads_vs_cnsd', 
    points = TRUE,
    create = function(){
        qcPlot(
            cnsdPlot,
            xcol = "cnsd",
            ycol = "total_num_reads",
            xlab = "sd(CN -> 1)",
            ylab = "Total # of Reads"
        )
        # abline(v = 0) # an example of direct plot manipulation
    }
)

# colData
# Classes ‘data.table’ and 'data.frame':  329 obs. of  38 variables:
#  $ barcode                    : chr [1:329(1d)] "AAACGGGTCTCGTGGG-1" "AAAGATGAGATCTGCT-1" "AACTCCCTCGATCCGG-1" "CCATGTCCAGCTATGT-1" ...
#  $ cell_id                    : chr  "0" "1" "10" "100" ...
#  $ effective_depth_of_coverage: num [1:329(1d)] 0.1249 0.0269 0.0586 0.094 0.047...
#  $ effective_reads_per_1Mbp   : int [1:329(1d)] 659 125 268 516 239 561 715 625 290 490 ...
#  $ est_cnv_resolution_mb      : num [1:329(1d)] 0.58 4.504 3.926 0.721 2.955 ...
#  $ frac_mapped_duplicates     : num [1:329(1d)] 0.139 0.126 0.129 0.14 0.121 ...
#  $ is_high_dimapd             : int [1:329(1d)] 0 0 1 0 0 0 0 0 0 0 ...
#  $ is_noisy                   : int [1:329(1d)] 0 0 1 0 1 0 0 0 0 0 ...
#  $ mean_ploidy                : num [1:329(1d)] 2.01 3.87 2.25 2.01 2.66 ...    
#  $ normalized_dimapd          : num [1:329(1d)] 1.3 1.5 2.08 1.26 1.7 ...       
#  $ normalized_mapd            : num [1:329(1d)] 0.0858 0.2547 0.241 0.0965 0.2055 ...
#  $ num_duplicate_reads        : int [1:329(1d)] 352062 68969 154580 266208 117157 282176 371025 325631 148296 248058 ...
#  $ num_lowmapq_reads          : int [1:329(1d)] 368158 136349 304821 224946 193821 329144 406299 372288 163224 272685 ...
#  $ num_mapped_dedup_reads     : int [1:329(1d)] 1797186 341058 729356 1405963 651362 1529771 1948322 1704397 790021 1334834 ...
#  $ num_unmapped_reads         : int [1:329(1d)] 8518 2210 5201 7857 4916 6641 9322 7418 5261 6659 ...
#  $ ploidy_confidence          : int [1:329(1d)] -2 -3 -4 -2 2 -2 -2 -2 8 -2 ... 
#  $ raw_dimapd                 : num [1:329(1d)] 1.3 1.5 2.08 1.26 1.7 ...       
#  $ raw_mapd                   : num [1:329(1d)] 0.0858 0.2547 0.241 0.0965 0.2055 ...
#  $ total_num_reads            : int [1:329(1d)] 2525924 548586 1193958 1904974 967256 2147732 2734968 2409734 1106802 1862236 ...
#------------------------------------------------------------------------
#  $ alignRate                  : num  0.851 0.747 0.74 0.878 0.795 ...
#  $ dupRate                    : num  0.164 0.168 0.175 0.159 0.152 ...
#  $ bad                        : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...      
#  $ keep                       : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ ploidy                     : int  2 2 2 2 2 2 2 2 2 2 ...
#  $ windowPower                : int  3 7 7 3 6 3 3 3 4 3 ...
#  $ modal_NA                   : int  2 2 2 2 2 2 2 2 2 2 ...
#  $ nrModelType                : chr  "singlePeak" "gain" "gain" "singlePeak" ...
#  $ q95_modal_NA               : num  74 187 410 60 189 64 81 71 70 57 ...       
#  $ halfWidthNA                : num  0.533 0.662 0.928 0.603 0.792 ...
#  $ ER_modal_NA                : num  58.4 140.5 280 46.1 135.4 ...
#  $ ER_ploidy                  : num  58.4 140.5 280 46.1 135.4 ...
#  $ RPA                        : num  29.2 70.3 140 23.1 67.7 ...
#  $ replicating                : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...   
#  $ fractionS                  : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ repGcRatio                 : num  NA 0.892 0.899 NA 0.892 ...
#  $ modelType                  : chr  "notReplicating" "notReplicating" "notReplicating" "notReplicating" ...
#  $ theta                      : num  367.96 15.01 7.88 297.13 11.14 ...

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
    # settings = settings$all_,
    # outcomes = reactive({ reactiveValuesToList(outcomes) }),
    # isReady  = reactive({ getStepReadiness(options$source, outcomes) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
