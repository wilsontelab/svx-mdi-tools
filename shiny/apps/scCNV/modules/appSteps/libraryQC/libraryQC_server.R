#----------------------------------------------------------------------
# static components to show plots with QC values and distributions
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
module <- 'libraryQC'   
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
outcomes <- reactiveValues() # outcomes[[sourceId]][i] <- TRUE if library i failed

#----------------------------------------------------------------------
# cascade from selected dataSource to its QC data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
project <- normalizeDataReactive(sourceId) # one Project has one or more Samples, each with one or more Cells

#----------------------------------------------------------------------
# make the QC plots
#----------------------------------------------------------------------
cellCol <- reactive({
    project <- project()
    req(project)
    d <- project$colData
    cols <- CONSTANTS$plotlyColors
    ifelse( d$bad,         cols$red, 
    ifelse(!d$keep,        cols$orange,
    ifelse( d$replicating, cols$blue,
                           cols$green)))
})
qcPlot2 <- function(xcol, ycol, xlab, ylab, limFn,
                    jitterX = FALSE, log10 = FALSE){
    project <- project()
    req(project)
    d <- project$colData
    x <- d[[xcol]]
    y <- d[[ycol]]
    if(log10) x <- log10(x)
    if(jitterX) x <- jitter(x)
    lim <- limFn(x, y)
    list(
        plotArgs = list(
            data.table(x = x, y = y),
            pch = 16,
            cex = 1,
            col = cellCol(),
            xlab = xlab,
            ylab = ylab,
            xaxs = "i",
            yaxs = "i"
        ),
        layout = list(
            width = 300,
            height = 250,
            pointsize = 9, 
            dpi = 96,
            mai = c(0.7, 0.7, 0.1, 0.1),
            xlim = lim$x, # OR can be read from plotArgs
            ylim = lim$y
        )
    )
}
#----------------------------------------------------------------------
totalReadLimits <- function(y) c(-100000, max(y, na.rm = TRUE) * 1.1)
fracPlotLimits <- function(x, y){
    list(
        x = c(min(x, na.rm = TRUE) - 0.05, max(x, na.rm = TRUE) + 0.05),
        y = totalReadLimits(y)
    )
}
wpPlotLimits <- function(x, y){
    list(
        x = c(min(x, na.rm = TRUE) - 0.5, max(x, na.rm = TRUE) + 0.5),
        y = totalReadLimits(y)
    )
}
cnsdPlotLimits <- function(x, y){
    list(
        x = c(min(x, na.rm = TRUE) - 0.05, max(x, na.rm = TRUE)* 1.1),
        y = totalReadLimits(y)
    )
}
#----------------------------------------------------------------------
arPlot <- mdiInteractivePlotServer(
    'nReads_vs_alignRate',      
    click = TRUE,
    contents = reactive({
        qcPlot2(
            xcol = "alignRate",
            ycol = "total_num_reads",
            xlab = "Alignment Rate",
            ylab = "Total # of Reads",
            limFn = fracPlotLimits
        )
    })
)
observeEvent(arPlot$click(), { setCellToPlot(
    arPlot$click(), 
    xcol = "alignRate",
    ycol = "total_num_reads",
    limFn = fracPlotLimits
)})
#----------------------------------------------------------------------
drPlot <- mdiInteractivePlotServer(
    'nReads_vs_dupRate',      
    click = TRUE,
    contents = reactive({
        qcPlot2(
            xcol = "dupRate",
            ycol = "total_num_reads",
            xlab = "Duplication Rate",
            ylab = "Total # of Reads",
            limFn = fracPlotLimits
        )
    })
)
observeEvent(drPlot$click(), { setCellToPlot(
    drPlot$click(), 
    xcol = "dupRate",
    ycol = "total_num_reads",
    limFn = fracPlotLimits
)})
#----------------------------------------------------------------------
wpPlot <- mdiInteractivePlotServer(
    'nReads_vs_windowPower',      
    click = FALSE,
    contents = reactive({
        qcPlot2(
            xcol = "windowPower",
            ycol = "total_num_reads",
            xlab = "Window Power",
            ylab = "Total # of Reads",
            limFn = wpPlotLimits, 
            jitterX = TRUE
        )
    })
)
#----------------------------------------------------------------------
cnsdPlot <- mdiInteractivePlotServer(
    'nReads_vs_cnsd',      
    click = TRUE,
    contents = reactive({
        qcPlot2(
            xcol = "cnsd",
            ycol = "total_num_reads",
            xlab = "sd(CN=1)",
            ylab = "Total # of Reads",
            limFn = cnsdPlotLimits
        )
    })
)
observeEvent(cnsdPlot$click(), { setCellToPlot(
    cnsdPlot$click(), 
    xcol = "cnsd",
    ycol = "total_num_reads",
    limFn = cnsdPlotLimits
)})

#----------------------------------------------------------------------
# plot one cell in response to QC plot click
#----------------------------------------------------------------------
setCellToPlot <- function(click, xcol, ycol, limFn){
    req(click)
    project <- project()
    req(project)
    d <- project$colData
    x <- d[[xcol]]
    y <- d[[ycol]]
    lim <- limFn(x, y)
    deltaX <- (click$coord$x - x) / diff(lim$x)
    deltaY <- (click$coord$y - y) / diff(lim$y)
    distance <- sqrt(deltaX ** 2 + deltaY ** 2)
    i <- which.min(distance)[1]
    cellToPlot(project$cells[[ d[i, cell_id] ]])
}
cellToPlot <- reactiveVal(NULL)
output$cellPlot <- renderUI({
    cell <- cellToPlot()
    req(cell)
    plotOneCellUI(project(), cell, settings)
})

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
