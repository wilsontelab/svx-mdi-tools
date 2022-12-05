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
outcomes <- reactiveValues()

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
    sourceData$colData
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

#----------------------------------------------------------------------
# parse genome windows for display, once per genome
#----------------------------------------------------------------------
windowPowers <- 0:7
minMappability <- 0.9 # TODO: expose as setting? if so, must add to cache key, below
parseWindows <- function(windowPower){
    windowSize <- as.integer(2 ** windowPower)
    windowSizeKey <- paste("w", windowSize, sep = "_")
    windows <- sourceData()$windows[[windowSizeKey]]
    i <- windows$mappability >= minMappability
    lastChromBins <- cumsum(windows[i, .N, by = "chrom"][[2]])
    list(
        i = i,
        N = sum(i),
        windows = windows[i],
        lastChromBins = lastChromBins,
        vLines = structure(
            c(0, lastChromBins) + 0.5, 
            color = rep('grey', length(lastChromBins) + 1)
        )
    )
}
mappableWindows <- reactive({
    sourceData <- sourceData()
    req(sourceData)    
    startSpinner(session, message = "parsing genome windows") 
    x <- windowsCache$get(
        "mappableWindows", 
        key = sourceData$metadata$assembly, 
        createFn = function(key, ...) lapply(windowPowers, parseWindows), 
        create = "always",
        from = "ram", 
        permanent = TRUE
    )   
    stopSpinner(session)
    x$value
})

#----------------------------------------------------------------------
# parse cell windows, once per window size
#----------------------------------------------------------------------
postProcessCellPower <- function(windowPower, ...){
    windowSize <- as.integer(2 ** windowPower)    
    mappableWindows <- mappableWindows()[[windowPower + 1]]
    cellData <- cellData()
    NR <- collapseVector(cellData$rawCounts, windowSize)[mappableWindows$i]
    NR_map <- NR / mappableWindows$windows$mappability
    ER_map_modal_CN <- peakValue(NR_map)
    list(
        windowPower = windowPower,
        windowSize = windowSize,
        minWindowPower = getMinWindowPower(cellData),
        NR_map = NR_map,
        ER_map_modal_CN = ER_map_modal_CN
        # TODO: initial fit? etc.
    )
}
cellWindowCounts <- reactive({
    req(input$windowPower)
    cellData <- cellData()
    req(cellData) 
    startSpinner(session, message = "post-processing cell") 
    x <- cellCache$get(
        "cellWindowCounts", 
        key = paste(cellData$colData$cell_id, "wp", input$windowPower, sep = "_"),
        createFn = postProcessCellPower, 
        create = "always",
        from = "disk", 
        permanent = TRUE,
        windowPower = input$windowPower
    )    
    stopSpinner(session)
    x$value
})
cellMappableWindows <- reactive({
    mappableWindows()[[cellWindowCounts()$windowPower + 1]]
})

# ----------------------------------------------------------------------
# handle cell inputs in top row
# ----------------------------------------------------------------------
observeEvent(input$windowPower, {
    req(input$windowPower)
    windowSize <- 2 ^ input$windowPower * CONSTANTS$binSize
    x <- if(windowSize > 1e6) paste(round(windowSize / 1e6, 1), "Mb") 
                         else paste(round(windowSize / 1e3, 1), "kb") 
    updateSliderInput(session, "windowPower", label = paste0("Window Size (", x, ")"))
})
observeEvent(cellWindowCounts(), {
    # dstr(cellWindowCounts())
#     dstr(cellData()$minWindowPower)
    updateSliderInput(session, "windowPower", min = cellWindowCounts()$minWindowPower)
})

# ----------------------------------------------------------------------
# GC bias plots
# ----------------------------------------------------------------------
NR_map_w_vs_gc_data <- reactive({
    cellWindowCounts <- cellWindowCounts()
    req(cellWindowCounts)
    windows <- cellMappableWindows()
    downsampleI <- sample.int(windows$N, min(2500, windows$N)) # TODO: expose downsample as setting
    data.table(
        x = windows$windows$gc_fraction[downsampleI],
        y = cellWindowCounts$NR_map[downsampleI]
    )
})
NR_map_w_vs_gc <- interactiveScatterplotServer(
    'NR_map_w_vs_gc',
    NR_map_w_vs_gc_data,
    accelerate = TRUE,
    xtitle = "Window Fraction GC",
    xrange = c(0.3, 0.6),
    # xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5)
    # selectable = 'h',
    # cacheReactive = preOutcomesSampleKey
)

#----------------------------------------------------------------------
# user node points set by clicking on upper genome plot
#----------------------------------------------------------------------
nodePointSize  <- 10
nodePointColor <- CONSTANTS$plotlyColors$orange
cnNodes <- data.table(x = integer(), y = integer())
observeEvent(NR_map_w_vs_bin$clicked(), {
    traceIndex <- 1  # 1 is the index of the overplot
    d <- NR_map_w_vs_bin$clicked()
    p <- plotlyProxy("NR_map_w_vs_bin-plotly", session, deferUntilFlush = FALSE)
    matchingNode <- cnNodes[, x == d$x]
    if(!any(matchingNode)){ # add a new node on first click
        cnNodes <<- rbind(cnNodes, data.table(
            x = d$x,
            y = d$y
        ))  
        plotlyProxyInvoke(p, "extendTraces", list(
            x = list(list(d$x)),
            y = list(list(d$y))
        ), list(traceIndex))            
    }
})

# ----------------------------------------------------------------------
# genome coordinate plots
# ----------------------------------------------------------------------
NR_map_w_vs_bin_data <- reactive({
    cellWindowCounts <- cellWindowCounts()
    req(cellWindowCounts)
    windows <- cellMappableWindows()
    data.table(
        x = 1:windows$N,
        y = cellWindowCounts$NR_map
    )
})
NR_map_w_vs_bin <- interactiveScatterplotServer(
    'NR_map_w_vs_bin',
    NR_map_w_vs_bin_data,
    overplot = data.table(x = -1000, y = 0),
    overplotMode = "markers",
    overplotColor = nodePointColor,
    overplotPointSize = nodePointSize,
    accelerate = TRUE,
    xtitle = "Genome Window",
    xrange = reactive({ c(1, cellMappableWindows()$N) }),
    xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5),
    hLines = function(...) NR_map_hLines(cellData, cellWindowCounts),
    vLines = reactive({ cellMappableWindows()$vLines }),
    clickable = TRUE
    # cacheReactive = preOutcomesSampleKey
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    # bm <- getModuleBookmark(id, module, bookmark, locks)
    # req(bm)
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
