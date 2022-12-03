#----------------------------------------------------------------------
# server components for the annotateCells appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
annotateCellsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'annotateCells'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# sample selection and loading
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
sourceData <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading source data")
    d <-loadSampleRds(sourceId)
    stopSpinner(session)
    d
})
colData <- reactive({ # i.e., the table of all cells
    sourceData <- sourceData()
    req(sourceData)
    # TODO: cache slopes
    slopes <- sapply(sourceData$colData$cell_id, function(cell_id){
        cell <- sourceData$cells[[cell_id]]
        if(is.null(cell$fit)) return(9999)
        gc <- cell$fit$gcFractions
        i <- gc >= 0.30 & gc <= 0.55
        gc <- gc[i]
        mu <- cell$fit$peak[i]
        coef(lm(mu ~ gc))[2]
    })
    sourceData$colData[order(-slopes)]
})
cellStepper <- listStepperButtonsServer( # inputs to pick a single sample or cell from the sample set
    id = 'cellStepper',
    dataReactive = colData,
    nameFn = function(x) x$cell_name
)
cellData <- reactive({
    sourceData <- sourceData()
    colData <- cellStepper$getCurrent()
    req(colData)
    list(
        colData = colData,
        rawCounts = sourceData$raw_counts[[colData$cell_id]]
    )
})
cellWindows <- reactive({
    windowSize <- as.integer(input$windowSize)
    windowSizeKey <- paste("w", windowSize, sep = "_")
    windows <- sourceData()$windows[[windowSizeKey]]
    i <- windows$mappability >= 0.9 # TODO: expose setting
    lastChromBins <- windows[i, .N, by = "chrom"][[2]]
    list(
        N = sum(i),
        windows = windows[i],
        NR_map = collapseVector(cellData()$rawCounts, windowSize)[i],
        lastChromBins = cumsum(lastChromBins)
    )
})

# ----------------------------------------------------------------------
# GC bias plots
# ----------------------------------------------------------------------
NR_map_w_vs_gc_data <- reactive({
    cellWindows <- cellWindows()
    i <- sample.int(cellWindows$N, min(2500, cellWindows$N)) # TODO: expose downsample as setting
    data.table(
        x = cellWindows$windows$gc_fraction[i],
        y = cellWindows$NR_map[i]
    )
})
NR_map_w_vs_gc <- interactiveScatterplotServer(
    'NR_map_w_vs_gc',
    NR_map_w_vs_gc_data,
    accelerate = TRUE,
    xtitle = "Window Fraction GC",
    xrange = c(0.3, 0.55),
    # xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5)
    # selectable = 'h',
    # cacheReactive = preOutcomesSampleKey
)

# ----------------------------------------------------------------------
# genome coordinate plots
# ----------------------------------------------------------------------
NR_map_w_vs_bin_data <- reactive({
    cellWindows <- cellWindows()
    data.table(
        x = 1:cellWindows$N,
        y = cellWindows$NR_map
    )
})
NR_map_w_vs_bin <- interactiveScatterplotServer(
    'NR_map_w_vs_bin',
    NR_map_w_vs_bin_data,
    accelerate = TRUE,
    xtitle = "Genome Window",
    # xrange = range_pad,
    # xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5),
    hLines = function(...) getPeakNR(cellWindows),
    vLines = function(...) getChromLines(cellWindows)
    # selectable = 'h',
    # cacheReactive = preOutcomesSampleKey
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
