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
        rawCounts = sourceData$raw_counts[[colData$cell_id]],
        cell = sourceData$cells[[colData$cell_id]]
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
    ploidy = 2,
    keep = TRUE
)
resetCellOutcomes <- function(){
    sourceId <- sourceId()
    cellData <- cellData()    
    req(sourceId, cellData)
    cell_id <- cellData$colData$cell_id
    defaults <- getCellDefaults(mappableWindows(), cellData)
    cellOutcomes[[sourceId]][[cell_id]] <- list(
        sourceId = sourceId,
        cell_id = cell_id,
        windowPower = defaults$windowPower,
        replicating = FALSE, # TODO: ? can this be inferred as default?
        modal_CN = cellData$colData$modal_CN,
        ploidy = 2, # TODO, inherit from pipeline, not sure if/where it was recorded
        keep = defaults$keep
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
        minWindowPower = getMinWindowPower(cellData), # used to prevent too-small choices of windowPower by user
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
        create = "always",
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
        setCellOutcome("windowPower", as.integer(input$windowPower))
        windowSize <- 2 ^ input$windowPower * CONSTANTS$binSize
        x <- if(windowSize > 1e6) paste(round(windowSize / 1e6, 1), "Mb") 
                             else paste(round(windowSize / 1e3, 1), "kb") 
        updateSliderInput(session, "windowPower", label = paste0("Window Size (", x, ")"))        
    }
})
observeEvent(input$replicating, {
    if(input$replicating != getCellOutcome("replicating")) {
        setCellOutcome("replicating", as.logical(input$replicating))
    } 
})
observeEvent(input$modal_CN, {
    if(input$modal_CN != getCellOutcome("modal_CN")) {
        setCellOutcome("modal_CN", as.integer(input$modal_CN))
    } 
})
observeEvent(input$ploidy, {
    if(input$ploidy != getCellOutcome("ploidy")) {
        setCellOutcome("ploidy", as.integer(input$ploidy))
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
    freezeReactiveValue(input, "ploidy")
    updateSliderInput(session, "ploidy", value = cellOutcomes$ploidy)
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
# user node points set by clicking on pre-fit genome plots
#----------------------------------------------------------------------
nodePointSize  <- 10
nullNodes <- list(
    cnv         = reactiveVal(data.table(x = -1000, y = 0)),
    replicating = reactiveVal(data.table(x = -1000, y = 0))
)
nodes <- nullNodes
observeEvent(NR_map_w_vs_bin0$clicked(), {
    # nodes$cnv <<- clickPlotNode(session, NR_map_w_vs_bin0, "NR_map_w_vs_bin0", nodes$cnv)
    nodes$cnv( clickPlotNode(session, NR_map_w_vs_bin0, "NR_map_w_vs_bin0", nodes$cnv()) )
})
observeEvent(byChromPlot$clicked(), {
    nodes$replicating( clickPlotNode(session, byChromPlot, "byChromPlot", nodes$replicating()) )
})
observeEvent(NR_map_w_bin_key0(), { # clear the nodes list whenever the corresponding plot changes
    nodes <<- nullNodes
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
            chrom = windows$windows$chrom,
            x = 1:windows$N,
            y = counts$NR_map
        ), #[chrom == "chrX"],
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
    overplot = nodes$cnv, #data.table(x = -1000, y = 0),
    overplotMode = "markers",
    overplotColor = CONSTANTS$plotlyColors$orange,
    overplotPointSize = nodePointSize,
    accelerate = TRUE,
    xtitle = "Genome Window",
    # xrange = function(d, ...) c(1, max(d$x)),
    xrange = function(d, ...) c(min(d$x), max(d$x)),
    xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5),
    yzeroline = FALSE,
    hLines = function(...) genomePlotBuffer0$hLines,
    vLines = function(...) genomePlotBuffer0$vLines,
    # vLines = function(...) 64970,
    clickable = TRUE,
    cacheReactive = NR_map_w_bin_key0
)

# ----------------------------------------------------------------------
# handle the by-chromosome panel to add replication nodes at fine level of granularity
# ----------------------------------------------------------------------
observeEvent(input$replicating, {
    toggle(
        selector = paste0("#", session$ns("byChromDiv")),
        condition = !is.null(input$replicating) & as.logical(input$replicating)
    )
})
chromStepper <- listStepperButtonsServer( # inputs to pick a single sample or cell from the sample set
    id = 'chromStepper',
    dataReactive = reactive({ cellMappableWindows()$windows[, .N, by = "chrom"] }),
    nameFn = function(x) x$chrom
)
NR_map_w_vs_bin_data_chrom <- reactive({
    chrom_ <- chromStepper$getCurrent()
    d <- NR_map_w_vs_bin_data0()    
    req(nrow(chrom_) > 0, d)
    d[chrom == chrom_$chrom]
})
byChromPlot <- interactiveScatterplotServer(
    'byChromPlot',
    NR_map_w_vs_bin_data_chrom,
    overplot = nodes$replicating,
    overplotMode = "markers",
    overplotColor = CONSTANTS$plotlyColors$red,
    overplotPointSize = nodePointSize,
    accelerate = TRUE,
    xtitle = "Genome Window",
    xrange = function(d, ...) range(d$x),
    xzeroline = FALSE,
    ytitle = '# Reads',
    yrange = function(...) range_pos(..., foldIQR = 2.5),
    yzeroline = FALSE,
    hLines = function(...) genomePlotBuffer0$hLines,
    clickable = TRUE
)

# ----------------------------------------------------------------------
# (re)fit a cells based on user outcomes
# ----------------------------------------------------------------------
output$density <- renderPlot({

    # dstr(colData())
    # dstr(cellData())
    # return(NULL)

    cell <- cellData()$cell
    dstr(cellData()$colData)

    return({
        plot(
            1:length(cell$cn), cell$cn, pch = ".", 
            ylim = c(0, 4) # c(0, cell$ER_modal_CN * 3)
        )
        lines(1:length(cell$cn), cell$hmm, col = "red")
    })



    # req(getCellOutcome("replicating"))
    modal_CN <- getCellOutcome("modal_CN")
    ploidy   <- getCellOutcome("ploidy")
    req(modal_CN == ploidy || modal_CN == 2 * ploidy)
    counts <- cellWindowCounts()
    ER_modal_CN <- counts$ER_map_modal_CN    
    NR <- counts$NR_map
    windowI <- !is.na(NR)
    NR <- NR[windowI]
    nBins <- length(NR)

    nReads <- sum(NR)

    peak_NR <- peakValue(NR)
    median_NR <- median(NR)

    # modalIsUnreplicated <- ER_modal_CN == ploidy
    # ER_replicated <- if(modalIsUnreplicated) 2 * ER_modal_CN else ER_modal_CN
    # ER_unreplicated <- ER_replicated / 2
    # sd_replicated <- (ER_replicated - ER_unreplicated) / (if(modalIsUnreplicated) 4 else 2)
    rpa <- peak_NR / ploidy

    windowChroms <- cellMappableWindows()$windows[windowI, chrom]
    chroms <- unique(windowChroms)
    chroms <- chroms[chroms != "chrY"]
    p <- sapply(chroms, function(chrom){
        median(NR[windowChroms == chrom])
        # wilcox.test(NR, NR[windowChroms == chrom])$p.value
    })

    cns <- 0:5
    getWorking <- function(v){
        q <- quantile(v, c(0.05, 0.95))
        v[v >= q[1] & v <= q[2] & v > 1]
    }
    NR_working <- getWorking(NR)
    rpa_working <- NR_working /  ploidy
    chromCn <- sapply(chroms, function(chrom){
        NR_chrom <- getWorking(NR[windowChroms == chrom])
        p <- sapply(cns, function(cn) wilcox.test(cn * rpa_working, NR_chrom, exact = FALSE)$p.value)
        cns[which.max(p)]
        # wilcox.test(NR, NR[windowChroms == chrom])$p.value
    })
    
    chromCn <- data.table(chrom = chroms, CN = chromCn)
    print(chromCn)

    plot(1:length(chroms), p, ylim = c(0, peak_NR + rpa * 2), 
         pch = 16, col = ifelse(p >= peak_NR - rpa / 2 & p <= peak_NR + rpa / 2, "blue", "red"))
    abline(h = peak_NR + rpa / 2 * -2:2)

    # agg <- aggregate(NR, list(NR), length)
    # agg$freq <- cumsum(agg[[2]]) / sum(agg[[2]])
    # plot(agg[[1]], agg$freq, typ = "l", ylim = c(0, 1))
    # for(chrom in chroms){
    #     chromNR <- NR[windowChroms == chrom]
    #     agg <- aggregate(chromNR, list(chromNR), length)
    #     agg$freq <- cumsum(agg[[2]]) / sum(agg[[2]])
    #     lines(agg[[1]], agg$freq)
    # }

    # modalIsUnreplicated <- ER_modal_CN == ploidy
    # ER_replicated <- if(modalIsUnreplicated) 2 * ER_modal_CN else ER_modal_CN
    # ER_unreplicated <- ER_replicated / 2
    # sd_replicated <- (ER_replicated - ER_unreplicated) / (if(modalIsUnreplicated) 4 else 2)
    # sd_unreplicated <- sd_replicated / 2
    # maxNR <- ER_replicated   + 2 * sd_replicated
    # minNR <- ER_unreplicated - 2 * sd_unreplicated
    # inRangeI <- NR >= minNR & NR <= maxNR
    # NR_inRange <- NR[inRangeI]
    # nInRange <- sum(inRangeI)

    # # eliminate outlier bins

    # # adjust the probability of being early replicating based on GC
    # d <- data.table(
    #     i = 1:nInRange,
    #     L_replicated   = dnorm(NR_inRange, ER_replicated,   sd = sd_replicated),
    #     L_unreplicated = dnorm(NR_inRange, ER_unreplicated, sd = sd_unreplicated) # here, i.e., dnorm() * pEarly_gc, etc. (very nearly just fracGC)
    # )

    # percentS <- 0:100
    # LL <- sapply(percentS / 100, function(fracS){
    #     d[, {
    #         x <- log(L_replicated * fracS + L_unreplicated * (1 - fracS))
    #         x[x == -Inf] <- -700 # suppress outliers            
    #         sum(x)
    #     }]
    # })
    # cellPercentS <- percentS[which.max(LL)]
    # nReplicatedBins <- round(NR_inRange * cellPercentS / 100, 0)

    # i <- d[order(-NR_inRange), i]
    # d[i, replicated := c(rep(TRUE, nReplicatedBins), rep(FALSE, nBins - nReplicatedBins))]

    # gc_fraction <- cellMappableWindows()$windows[windowI, gc_fraction][inRangeI]
    # dprint(median(gc_fraction[d$replicated == TRUE]))
    # dprint(median(gc_fraction[d$replicated == FALSE]))

    # plot(1:nInRange, NR_inRange, pch = 16, cex = 0.5, col = ifelse(d$replicated, "red", "blue"))
    # abline(h = ER_replicated   + 0 * sd_replicated, col = "black")
    # abline(h = ER_unreplicated - 0 * sd_unreplicated, col = "black")
    # abline(h = ER_replicated   + 2 * sd_replicated, col = "red")
    # abline(h = ER_unreplicated - 2 * sd_unreplicated, col = "red")
    # plot(1:101, LL, ylim = c(quantile(LL, 0.1), max(LL)))
    # plot(density(NR))
})

    # CELLS
    # aneuploid = 91 (gain), 81 (loss), 76 (massive amplification, CN=165)
    # replicating = 24 (late, clean), 32 (early, XY), 82 (late, with aneuploidy), 79 (early), 
    #               57 (75% clean), 68 (late), 69 (mid, XY, CNVs), 79? (early?), 84 (mid, with aneuploidy, bad QC level)
    # segments = 77 (multiple)

observeEvent(input$fitCell, {

    dmsg("cell fit pending")
    dstr(nodes)

    counts <- cellWindowCounts()
    modal_CN <- getCellOutcome("modal_CN")
    working <- data.table(
        chrom = cellMappableWindows()$chrom,
        NR_map = counts$NR_map,
        node = FALSE
    )
    working$node[nodes$x] <- TRUE

    dstr(working)


    
    # https://stephens999.github.io/fiveMinuteStats/intro_to_em.html

    # TODO: handle chrX specially, based on median NR_map

    # # use nodes to set CN per chromosome per segment    
    # working[, cn := { 
    #     if(sum(node) == 0) modal_CN # chroms with no nodes set as modal CN
    #     else if(sum(node) == 1) round(NR_map[node] / counts$ER_map_modal_CN * modal_CN, 0)
    #     else {
    #         nodeCn <- round(NR_map[node] / counts$ER_map_modal_CN * modal_CN, 0)
    #         segmentEnds <- sapply(1:(sum(node) - 1), function(i){
    #             node1 <- which(node[i])
    #             node2 <- which(node[i + 1])
    #             NR <- NR_map[node1:node2]
    #             delta <- sapply(node1:(node2 - 1), function(j){
    #                 left <- NR_map[node1:j] - nodeCn[i]
    #                 right <- NR_map[(j + 1):node2] - nodeCn[i + 1]
    #                 sum(left + right)
    #             })
    #             which(delta == min(delta))[1]

    #         })

    #     }
    # # use chrom-level hmm (or other?) to find the single breakpoint between segments    
    # }, by = chrom]

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
