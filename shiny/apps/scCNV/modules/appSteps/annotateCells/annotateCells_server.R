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
sourceId <- dataSourceTableServer("source", selection = "single")
sourceData <- reactive({
    toggle(selector = "#annotateCellsInputs", condition = FALSE)
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading source data")    
    d <-loadSampleRds(sourceId)    
    toggle(selector = "#annotateCellsInputs", condition = TRUE)
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
# cell-level module outcomes, created by user clicks
#----------------------------------------------------------------------
cellOutcomes <- reactiveValues() # the primary output of this module = user annotation of cells
nullCellOutcomes <- list(
    sourceId = NA,
    cell_id = NA,
    windowPower = 2,
    replicating = FALSE,
    modal_CN = 2,
    keep = TRUE
)
resetCellOutcomes <- function(){
    sourceId <- sourceId()
    cellData <- cellData()    
    req(sourceId, cellData)
    cell_id <- cellData$colData$cell_id
    cellOutcomes[[sourceId]][[cell_id]] <- list(
        sourceId = sourceId,
        cell_id = cell_id,
        windowPower = getMinWindowPower(cellData),
        replicating = FALSE,
        modal_CN = cellData$colData$modal_CN,
        keep = TRUE
    )
}
getCellOutcomes <- function(sourceId = NULL, cell_id = NULL){
    if(is.null(sourceId)) sourceId <- sourceId()
    if(is.null(cell_id)) cell_id <- cellData()$colData$cell_id  
    req(sourceId, cell_id)
    if(is.null(cellOutcomes[[sourceId]][[cell_id]])) resetCellOutcomes()
    cellOutcomes[[sourceId]][[cell_id]]
}
getCellOutcome <- function(key, reactive = TRUE){
    if(reactive) getCellOutcomes()[[key]]
    else isolate({ getCellOutcomes()[[key]] })
}
setCellOutcome <- function(key, value){
    sourceId <- sourceId()
    cellData <- cellData()    
    req(sourceId, cellData)
    cell_id <- cellData$colData$cell_id
    cO <- getCellOutcomes(sourceId, cell_id)
    cO[[key]] <- value
    cellOutcomes[[sourceId]][[cell_id]] <- cO
}
keptCells <- reactive({
    NULL
})

#----------------------------------------------------------------------
# parse genome windows for display, once per genome
#----------------------------------------------------------------------
parseWindows <- function(windowPower, windows){
    windowSize <- as.integer(2 ** windowPower)
    windowSizeKey <- paste("w", windowSize, sep = "_")
    windows <- windows[[windowSizeKey]]
    i <- windows$mappability >= CONSTANTS$minMappability
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
        createFn = function(key, ...) lapply(CONSTANTS$windowPowers, parseWindows, sourceData$windows), 
        create = "asNeeded",
        from = "ram", 
        permanent = TRUE
    )   
    stopSpinner(session)
    x$value
})

#----------------------------------------------------------------------
# parse cell windows, once per window size
#----------------------------------------------------------------------
postProcessCellPower <- function(cellData, mappableWindows, windowPower, ...){
    mappableWindows <- mappableWindows[[windowPower + 1]]
    windowSize <- as.integer(2 ** windowPower)    
    NR <- collapseVector(cellData$rawCounts, windowSize)[mappableWindows$i]
    NR_map <- NR / mappableWindows$windows$mappability
    ER_map_modal_CN <- peakValue(NR_map)
    list(
        colData = cellData$colData,
        windowPower = windowPower,
        windowSize = windowSize,
        minWindowPower = getMinWindowPower(cellData),
        NR_map = NR_map,
        ER_map_modal_CN = ER_map_modal_CN
        # TODO: initial fit? etc.
    )
}
cellWindowCounts <- reactive({
    cellData <- cellData()
    mappableWindows <- mappableWindows()
    req(cellData, mappableWindows) 
    windowPower <- getCellOutcome("windowPower")
    startSpinner(session, message = "post-processing cell") 
    x <- cellCache$get(
        "cellWindowCounts", 
        key = paste(cellData$colData$cell_id, "wp", windowPower, sep = "_"),
        createFn = postProcessCellPower, 
        create = "asNeeded",
        from = "disk", 
        permanent = TRUE,
        cellData = cellData,
        mappableWindows = mappableWindows,
        windowPower = windowPower
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
    minWindowPower <- cellWindowCounts()$minWindowPower
    if(input$windowPower < minWindowPower){
        updateSliderInput(session, "windowPower", value = minWindowPower)
    } else if(input$windowPower != getCellOutcome("windowPower")) {
        setCellOutcome("windowPower", input$windowPower)
        windowSize <- 2 ^ input$windowPower * CONSTANTS$binSize
        x <- if(windowSize > 1e6) paste(round(windowSize / 1e6, 1), "Mb") 
                             else paste(round(windowSize / 1e3, 1), "kb") 
        updateSliderInput(session, "windowPower", label = paste0("Window Size (", x, ")"))        
    }
})
observeEvent(input$replicating, {
    if(input$replicating != getCellOutcome("replicating")) {
        setCellOutcome("replicating", input$replicating)
    } 
})
observeEvent(input$modal_CN, {
    if(input$modal_CN != getCellOutcome("modal_CN")) {
        setCellOutcome("modal_CN", input$modal_CN)
    } 
})
updateKeepRejectCellButton <- function(keep){
    label <- if(keep) "Reject" else "Keep"
    style <- if(keep) "warning" else "success"
    updateButton(session, session$ns("keepRejectCell"), label = label, style = style)    
}
observeEvent(input$keepRejectCell, {
    keep <- !getCellOutcomes()$keep
    setCellOutcome("keep", keep)
    updateKeepRejectCellButton(keep)
})
observeEvent(input$resetCell, {
    cO <- getCellOutcomes()
    match <- sapply(names(nullCellOutcomes), function(x) cO[[x]] == nullCellOutcomes[[x]])
    if(all(match)) invalidateNodePlot( invalidateNodePlot() + 1 ) # not working, likely due to plot cache
    else updateCellInputs(nullCellOutcomes)
})
updateCellInputs <- function(cellOutcomes){
    freezeReactiveValue(input, "windowPower")
    updateSliderInput(session, "windowPower", value = cellOutcomes$windowPower)        
    freezeReactiveValue(input, "replicating")
    updateRadioButtons(session, "replicating", selected = cellOutcomes$replicating)
    freezeReactiveValue(input, "modal_CN")
    updateSliderInput(session, "modal_CN", value = cellOutcomes$modal_CN)
    updateKeepRejectCellButton(cellOutcomes$keep)     
}
observeEvent(getCellOutcomes(), {
    cO <- getCellOutcomes()
    req(cO)
    updateCellInputs(cO)
})

# ----------------------------------------------------------------------
# define keys for triggering and caching interactive plots
# ----------------------------------------------------------------------
NR_map_w_gc_key0 <- reactive({    
    getCellOutcomes()[c('sourceId', 'cell_id', 'windowPower')]
})
NR_map_w_bin_key0 <- reactive({ 
    getCellOutcomes()[c('sourceId', 'cell_id', 'windowPower', 'modal_CN')]
})

# ----------------------------------------------------------------------
# GC bias plot, pre-fit
# ----------------------------------------------------------------------
invalidateGcPlot0 <- reactiveVal(0)
gcPlotBuffer0 <- list()
observeEvent({
    NR_map_w_gc_key0()
}, {
    windows <- cellMappableWindows()    
    cellWindowCounts <- cellWindowCounts()
    downsampleI <- sample.int(windows$N, min(2500, windows$N)) # TODO: expose downsample as setting
    gcPlotBuffer0 <<- list(
        points = data.table(
            x = windows$windows$gc_fraction[downsampleI],
            y = cellWindowCounts$NR_map[downsampleI]
        )
    )
    invalidateGcPlot0( invalidateGcPlot0() + 1 )
})
NR_map_w_vs_gc_data0 <- reactive({
    invalidateGcPlot0()
    gcPlotBuffer0$points
})
NR_map_w_vs_gc0 <- interactiveScatterplotServer(
    'NR_map_w_vs_gc0',
    NR_map_w_vs_gc_data0,
    accelerate = TRUE,
    xtitle = "Window Fraction GC",
    xrange = c(0.3, 0.6),
    xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5),
    yzeroline = FALSE,
    cacheReactive = NR_map_w_gc_key0
)

#----------------------------------------------------------------------
# user node points set by clicking on upper genome plot
#----------------------------------------------------------------------
nodePointSize  <- 10
nodePointColor <- CONSTANTS$plotlyColors$orange
nullNodes <- data.table(x = integer(), y = integer())
pendingNodes <- nullNodes
addNode <- Vectorize(function(x, y){
    p <- plotlyProxy("NR_map_w_vs_bin0-plotly", session, deferUntilFlush = FALSE)      
    plotlyProxyInvoke(p, "extendTraces", list(
        x = list(list(x)),
        y = list(list(y))
    ), list(1)) # 1 is the index of the overplot trace 
})
observeEvent(NR_map_w_vs_bin0$clicked(), {
    d <- NR_map_w_vs_bin0$clicked()
    matchingNode <- pendingNodes[, x == d$x]
    if(!any(matchingNode)) {
        addNode(d$x, d$y)
        pendingNodes <<- rbind(pendingNodes, data.table(x = d$x, y = d$y)) 
    }
})
observeEvent(NR_map_w_bin_key0(), { # clear the nodes list whenever the corresponding plot changes
    pendingNodes <<- nullNodes
})

# ----------------------------------------------------------------------
# genome coordinate plots
# ----------------------------------------------------------------------
invalidateGenomePlot0 <- reactiveVal(0)
genomePlotBuffer0 <- list()
observeEvent({
    NR_map_w_bin_key0()
}, {
    windows <- cellMappableWindows()
    counts <- cellWindowCounts()
    cO <- getCellOutcomes()
    genomePlotBuffer0 <<- list(
        points = data.table(
            x = 1:windows$N,
            y = counts$NR_map
        ),
        hLines = NR_map_hLines(cO$modal_CN, counts),
        vLines = windows$vLines
    )
    invalidateGenomePlot0( invalidateGenomePlot0() + 1 )
})
NR_map_w_vs_bin_data0 <- reactive({
    invalidateGenomePlot0()
    genomePlotBuffer0$points
})
NR_map_w_vs_bin0 <- interactiveScatterplotServer(
    'NR_map_w_vs_bin0',
    NR_map_w_vs_bin_data0,
    overplot = data.table(x = -1000, y = 0),
    overplotMode = "markers",
    overplotColor = nodePointColor,
    overplotPointSize = nodePointSize,
    accelerate = TRUE,
    xtitle = "Genome Window",
    xrange = function(d, ...) c(1, max(d$x)),
    xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5),
    yzeroline = FALSE,
    hLines = function(...) genomePlotBuffer0$hLines,
    vLines = function(...) genomePlotBuffer0$vLines,
    clickable = TRUE,
    cacheReactive = NR_map_w_bin_key0
)

# ----------------------------------------------------------------------
# (re)fit a cells based on user outcomes
# ----------------------------------------------------------------------
observeEvent(input$fitCell, {

    dmsg("cell fit pending")
    dstr(pendingNodes)
    # use nodes to set CN per chromosome per segment
    # use chrom-level hmm (or other?) to find the single breakpoint between segments
    # chroms with no nodes are assumed to be modal_CN
    # chroms with one node are assumed to be the nearest multiple to node Y

    # for replicating cells, solve a mixture model for CN = 2 and CN = 4

    # use both of the above to adjust the nb gc curve fit

    # use the nb gc curve fit to replot the cell
    # and solve hmm

})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    cO <- bm$outcomes$cellOutcomes
    if(!is.null(cO)) cellOutcomes <<- listToReactiveValues(cO)
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(
        cellOutcomes = reactive({ reactiveValuesToList(cellOutcomes) })
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }), # TODO: all cells fitted?
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
