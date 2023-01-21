#----------------------------------------------------------------------
# server components for the keepRejectCells appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
keepRejectCellsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'keepRejectCells'   
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
overrides <- reactiveValues()
isUserOverride <- Vectorize(function(cell_id, sourceId = NULL) {
    if(is.null(sourceId)) sourceId <- sourceId()
    !is.null(overrides[[sourceId]][[cell_id]]) && overrides[[sourceId]][[cell_id]]
})
getKeep <- function(bad, keep, cell_id, sourceId) {
    overridden <- isUserOverride(cell_id, sourceId)
    (!bad &  keep & !overridden) | 
    (!bad & !keep &  overridden)
}
getRejected <- function(bad, keep, cell_id, sourceId){
    overridden <- isUserOverride(cell_id, sourceId)
    (!bad & !keep & !overridden) | 
    (!bad &  keep &  overridden)
}

#----------------------------------------------------------------------
# source data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
# projectName <- projectNameReactive(sourceId)
sample <- normalizeDataReactive(sourceId)
manifest <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    fread(getSourceFilePath(sourceId, "manifestFile"))
})
cells <- reactive({
    sample <- sample()
    req(sample)  
    sourceId <- sourceId()
    sampleCellIds <- sampleCellIds()
    i <- sample$colData[, {
        switch(
            input$cellStatus,
            All    = rep(TRUE, .N),
            Bad    = bad,
            Keep   = getKeep(bad, keep, cell_id, sourceId),
            Reject = getRejected(bad, keep, cell_id, sourceId)
        ) & switch(
            input$replicating,
            All = TRUE,
            Yes = !bad &  replicating,
            No  = !bad & !replicating
        ) & (sampleCellIds[1] == ALL | cell_id %in% sampleCellIds) &
          if(input$cellIdFilter != ALL) cell_id == input$cellIdFilter else TRUE
    }]
    n <- sum(i)
    list(
        i = i,
        n = n,
        maxPage = ceiling(n / input$cellsPerPage)
    )
})
slopes <- reactive({
    sample <- sample()
    req(sample)    
    sapply(sample$colData$cell_id, function(cell_id){
        cell <- sample$cells[[cell_id]]
        if(cell$badCell) return(NA)
        shapeModel <- settings$get("Page_Options", "Shape_Model")
        cw <- cell$windows[[cell$windowPower + 1]][[shapeModel]]

        dstr(cw)
        gc_fit <- cw$sequential$gc_fit
        if(is.null(gc_fit)) return(NA)
        gc <- gc_fit$gcFractions
        i <- gc >= 0.35 & gc <= 0.5
        gc <- gc[i]
        mu <- gc_fit$peak[i] / cw$RPA
        coef(lm(mu ~ gc))[2]
    })
})

#----------------------------------------------------------------------
# cascade update sample and cell selectors
#----------------------------------------------------------------------
ALL <- "__ALL__"
sampleCellIds <- reactiveVal(ALL)
updateCellFilter <- function(){
    manifest <- manifest()
    req(manifest)
    i <- manifest[, input$sampleNameFilter == ALL | Sample_Name == input$sampleNameFilter]
    choices <- manifest[i, Sample_ID]
    names(choices) <- manifest[i, Description]
    sampleCellIds(choices)
    freezeReactiveValue(input, "pageNumber")
    updateTextInput(session, "pageNumber", value = 1)
    freezeReactiveValue(input, "cellIdFilter")
    updateSelectInput(session, "cellIdFilter", choices = c(All = ALL, sort(choices)), selected = ALL)    
}
observeEvent(sourceId(), {
    manifest <- manifest()
    req(manifest)
    choices <- manifest$Sample_Name
    names(choices) <- choices
    freezeReactiveValue(input, "sampleNameFilter")
    updateSelectInput(session, "sampleNameFilter", choices = c(All = ALL, sort(unique(choices))), selected = ALL)
    updateCellFilter()
})
observeEvent(input$sampleNameFilter, {
    updateCellFilter()
})
observeEvent(cells(), {
    freezeReactiveValue(input, "pageNumber")
    updateTextInput(session, "pageNumber", value = 1)
})

#----------------------------------------------------------------------
# a stack of individual cell plots
#----------------------------------------------------------------------
pointOpacity <- function(x){
    min(1, max(0.15, -1/6000 * length(x) + 1))
}
getCnNAColor <- function(cnNa, pointOpacity = NULL){ # by CN (regardless of ploidy)
    if(is.null(pointOpacity)) pointOpacity <- pointOpacity(cnNa)
    defaultPointColor <- rgb(0, 0, 0, pointOpacity)
    cnNaCols <- c(     
        defaultPointColor,                # 0 = black/grey (absence of color/copies...)
        rgb(0,   0,    1,   pointOpacity),       # 1 = blue ("cool" colors are losses)
        rgb(0.1, 0.8,  0.1, pointOpacity), # 2 = subdued green ("good/go" for typical CN neutral)
        rgb(1,   0,    0,   pointOpacity),       # 3 = red ("hot" colors are gains)
        rgb(1,   0.65, 0,   pointOpacity),    # 4 = orange
        defaultPointColor                 # 5 = back to black/grey to make it obvious
    )
    col <- cnNaCols[cnNa + 1]
    col[is.na(col)] <- defaultPointColor
    col
}
plotCellByWindow <- function(y, ylab, h, v, col = NULL, yaxt = NULL, ymax = NULL){
    if(is.null(col)) col <- rgb(0, 0, 0, pointOpacity(y))
    x <- 1:length(y)
    if(is.null(ymax)) ymax <- max(h)
    plot(NA, NA, bty = "n", xaxt = "n", yaxt = yaxt,
         xlab = NULL, ylab = ylab, xlim = range(x), ylim = c(0, ymax))
    abline(h = h, v = v, col = "grey")
    points(x, y, pch = 19, cex = 0.4, col = col)
}
plotCellQC <- function(pngFile, sample, cell){
    png(
        filename = pngFile,
        width  = 1.5 * 6, 
        height = 1.35, 
        units = "in", 
        pointsize = 7,
        bg = "white",
        res = 96, # i.e., optimized for on-screen display
        type = "cairo"
    )
    layout(matrix(c(c(1,1), rep(c(2,3), 5)), nrow = 2, ncol = 6))

    w <- sample$windows[[cell$windowPower + 1]]
    v <- c(0, w[, max(i), by = "chrom"][[2]]) + 0.5

    if(cell$badCell){
        cw <- cell$windows$unshaped # bad cells only have a windows entry for windowPower 7, which is unshaped

        # plot NR_wms vs. gc_w
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        N <- length(cw$NR_wms)
        I <- sample.int(N, min(1e4, N))
        plot(w$gc_fraction[I], cw$NR_wms[I], 
            xlab = "Fraction GC", xlim = c(0.3, 0.6), 
            ylab = "# of Reads", #ylim = c(0, quantile(RPA, 0.975, na.rm = TRUE) * 1.5),
            pch = 19, cex = 0.4)

        # plot NR_wm vs. window index, i.e., pre-normalization input
        par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
        plotCellByWindow(cw$NR_wms, "# Reads", h = max(cw$NR_wms, na.rm = TRUE), v = v)
        
    } else {
        shapeModel <- settings$get("Page_Options", "Shape_Model")
        replicationModel <- settings$get("Page_Options", "Replication_Model")
        forceSequential <- cell$cellIsReplicating && replicationModel == "Sequential"
        shapeKey <- if(shapeModel == "Unshaped") "unshaped" else "shaped"
        repKey <- if(!cell$cellIsReplicating || forceSequential) "sequential" else "composite"
        cw <- cell$windows[[cell$windowPower + 1]][[shapeKey]]
        cww <- cw[[repKey]]
        cww$NA_ <- cww$HMM + cww$NAR

        # plot NR_wms vs. gc_w
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        N <- length(cw$NR_wms)
        I <- sample.int(N, min(1e4, N))
        RPA <- cw$NR_wms[I] / (if(forceSequential) cww$HMM[I] else cww$NA_[I])
        plot(w$gc_fraction[I], RPA, 
            xlab = "Fraction GC", xlim = c(0.3, 0.6), 
            ylab = "Reads Per Allele", ylim = c(0, quantile(RPA, 0.975, na.rm = TRUE) * 1.5),
            pch = 19, cex = 0.4, col =  if(forceSequential) rgb(0, 0, 0, 0.05) else getCnNAColor(cww$NA_[I], 0.05))
        gc_fit <- if(forceSequential) cww$gc_fit else cell$replicationModel[[shapeKey]]$gc_fit 
        with(gc_fit, { for(percentile in c(0.025, 0.975)) lines(
            gcFractions, 
            qnbinom(percentile, size = theta, mu = mu), 
            lty = 3, lwd = 1.5, col = "grey30"
        ) })
        
        # plot NR_wms vs. window index, i.e., pre-normalization input
        par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
        maxPlotNa <- max(cww$NA_, na.rm = TRUE) + 1.4
        plotCellByWindow(cw$NR_wms, "# Reads", h = cw$RPA * 0:round(maxPlotNa, 0), v = v, ymax = cw$RPA * maxPlotNa)

        # plot CN vs. window index, i.e., post-normalization
        par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
        maxPlotNa <- max(cww$HMM, na.rm = TRUE) + 1
        plotCellByWindow(cww$CN, "CN", h = 0:maxPlotNa, v = v, col = getCnNAColor(cww$HMM), yaxt = "n")
        lines(1:length(cww$HMM), cww$HMM, col = "red")  
        axis(2, at=0:maxPlotNa, labels=0:maxPlotNa)      
    }

    dev.off()
}
getCellCompositePlot <- function(sample, cell){
    fileName <- paste(cell$cell_id, "qc", "png", sep = ".")
    pngFile <- file.path(sample$qcPlotsDir, fileName)

    if(TRUE || !file.exists(pngFile)) plotCellQC(pngFile, sample, cell)

    tags$img(src = pngFileToBase64(pngFile), style = "vertical-align: top;")
}
getCellSummary <- function(colData, sourceId){
    manifest <- manifest()
    desc <- manifest[as.character(Sample_ID) == colData$cell_id, Description]
    prefix <- session$ns("")
    kept <- colData[, getKeep(bad, keep, cell_id, sourceId)]
    override <- if(colData$bad) "" else if(colData$keep){
        if( kept) TRUE else FALSE        
    } else {
        if(!kept) TRUE else FALSE
    }
    tags$div(
        style = "display: inline-block; padding: 5px; margin-left:5px;",
        tags$div(tags$strong( paste("cell:", colData$cell_id, '(', desc, ')') )), 
        tags$div(paste("windowPower: ", colData$windowPower, '(', commify(2**colData$windowPower * 20), 'kb )')),
        if(colData$bad) tags$div("bad cell, not processed") else tagList(
            tags$div(paste("repGcRatio:", round(colData$repGcRatio, 2), '(', colData$modelType, ',', round(colData$fractionS, 2), ')')), 
            tags$div(paste("cnsd:", round(colData$cnsd, 2), '(', if(colData$keep) "" else "not", "passed", ')')),
            tags$div(
                style = "margin-top: 5px;",
                tags$button(
                    onclick = paste0("cellKeepReject('", prefix, "', '", colData$cell_id, "', '", override, "')"),
                    if(kept) "Reject" else "Keep"
                )   
            )        
        ) 
    )
}
output$cellPlots <- renderUI({ 
    sourceId <- sourceId()     
    cells <- cells()
    req(cells)
    if(cells$n == 0) return("no cells to plot")
    startSpinner(session, message = "loading cell plots")
    sample <- sample()
    sortBy <- settings$get("Page_Options", "Sort_By")   
    ascDesc <- if(settings$get("Page_Options", "Order") == "Ascending") 1 else -1
    colData <- sample$colData[cells$i][order(ascDesc * switch(
        sortBy,
        windowPower = sample$colData[cells$i, windowPower + bad],
        # manifest = manifest()[Sample_ID %in% sample$colData[cells$i][, cell_id], .I],
        cnsd = sample$colData[cells$i, cnsd],
        fractionS = sample$colData[cells$i, fractionS],
        slope = slopes()[cells$i],
        sample$colData[cells$i, cell_id]
    ))]
    cellIndices <- 1:input$cellsPerPage + input$cellsPerPage * (as.integer(input$pageNumber) - 1)
    x <- lapply(cellIndices, function(i){
        if(i > cells$n) return(NULL)
        colData <- colData[i, ]
        cell <- sample$cells[[colData$cell_id]]
        tags$div(
            style = paste("border: 1px solid grey; white-space: nowrap;"),
            tagList(
                getCellCompositePlot(sample, cell),
                getCellSummary(colData, sourceId)                
            )
        )
    })
    stopSpinner(session)
    x
})

#----------------------------------------------------------------------
# pagination control
#----------------------------------------------------------------------
output$nCellsInType <- renderText({
    cells <- cells()
    req(cells)
    paste0(cells$n, " cells (", cells$maxPage, " pages)")
})
observeEvent(input$prevPage, {
    pageNumber <- as.integer(input$pageNumber)
    if(pageNumber == 1) return()
    updateTextInput(session, "pageNumber", value = pageNumber - 1)
})
observeEvent(input$nextPage, {
    pageNumber <- as.integer(input$pageNumber)
    if(pageNumber == cells()$maxPage) return()
    updateTextInput(session, "pageNumber", value = pageNumber + 1)
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    bo <- bm$outcomes

    # for(sourceId in names(bo$userOverrides)) userOverrides[[sourceId]] <- bo$userOverrides[[sourceId]]
    # for(sourceId in names(bo$userModalCN))   userModalCN[[sourceId]]   <- bo$userModalCN[[sourceId]]

    settings$replace(bm$settings)
    # updateTextInput(session, 'cellsPerPage', value = bm$input$cellsPerPage)
    # updateTextInput(session, 'pageNumber',   value = bm$input$pageNumber)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,    
    NULL
    # outcomes = list(
    #     userOverrides = userOverrides,
    #     userModalCN   = userModalCN
    # ),
    # # isReady = reactive({ getStepReadiness(options$source, ...) }),
    # getSampleFilePrefix = function(sourceId, keep = NULL, colData = NULL){
    #     if(is.null(keep)){
    #         if(is.null(colData)) colData <- {
    #             dataFilePath <- getSourceFilePath(sourceId, "normalizeFile")
    #             readRDS(dataFilePath)$colData
    #         }
    #         keep <- getKeep(colData$rejected, colData$cell_id, sourceId)            
    #     }
    #     key <- digest(keep)
    #     expandSourceFilePath(sourceId, key)
    # },
    # sampleSummaryColumns = reactive({
    #     sources <- getStepReturnValueByType('upload', 'outcomes')$sources()
    #     startSpinner(session, message = "get sample summaries") 
    #     x <- do.call(rbind, lapply(names(sources), function(sourceId){
    #         normalizeFilePath <- getSourceFilePath(sourceId, "normalizeFile")
    #         colData <- readRDS(normalizeFilePath)$colData
    #         keep <- getKeep(colData$rejected, colData$cell_id, sourceId)
    #         batchFilePath <- paste(app$adjust$getSampleFilePrefix(sourceId, keep), "batchNormalized.rds", sep = ".")
    #         batchFileExists <- file.exists(batchFilePath)
    #         colData[, .(    
    #             sourceId = sourceId,
    #             normalizeFilePath = normalizeFilePath,
    #             batchFilePath = batchFilePath,
    #             batchFileExists = batchFileExists,
    #             cell_id = list(colData$cell_id[keep]),
    #             cells = .N,
    #             keep     = sum(keep),
    #             reject   = sum(getRejected(rejected, cell_id, sourceId)),
    #             override = sum(isUserOverride(cell_id, sourceId)),
    #             status = if(batchFileExists) 
    #                 '<b style="margin-left: 10px;"><i class="fa-solid fa-check fa-lg" style="color: #0a0;"></i><b>' 
    #                 else "pending"
    #         )]   
    #     }))
    #     stopSpinner(session)
    #     x
    # })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
