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
# initialize module and outcomes
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
    size = "m"
)
outcomes <- reactiveValues() # outcomes[[sourceId]][i] <- TRUE if library i failed
overrides <- reactiveValues()
keptCnvs <- reactiveValues()
showFiltersAsOverrides <- reactive({
    x <- settings$get("Page_Options", "Show_Filters_As")
    if(isTruthy(x)) x == "User Overrides" else TRUE
})
isUserOverride <- Vectorize(function(cell_id, key, sourceId, forceOverrides = TRUE) {
    if(!forceOverrides && !showFiltersAsOverrides()) FALSE 
    else !is.null(overrides[[sourceId]][[cell_id]][[key]]) && overrides[[sourceId]][[cell_id]][[key]]
})
getKeep <- function(bad, keep, cell_id, sourceId = NULL, forceOverrides = TRUE) {
    if(is.null(sourceId)) sourceId <- sourceId()
    overridden <- isUserOverride(cell_id, "keep", sourceId, forceOverrides)
    (!bad &  keep & !overridden) | 
    (!bad & !keep &  overridden)
}
getRejected <- function(bad, keep, cell_id, sourceId = NULL, forceOverrides = TRUE){
    if(is.null(sourceId)) sourceId <- sourceId()
    overridden <- isUserOverride(cell_id, "keep", sourceId, forceOverrides)
    (!bad & !keep & !overridden) | 
    (!bad &  keep &  overridden)
}
getReplicating <- function(bad, replicating, cell_id, sourceId = NULL, forceOverrides = TRUE){
    if(is.null(sourceId)) sourceId <- sourceId()
    overridden <- isUserOverride(cell_id, "replicating", sourceId, forceOverrides)
    !bad & replicating & !overridden
}
toggleKeptCnv <- function(cnv, sourceId = NULL){
    if(is.null(sourceId)) sourceId <- sourceId()
    key <- getCnvKey(cnv, sourceId)
    keptCnvs[[key]] <- if(is.null(keptCnvs[[key]])) TRUE else NULL
}
isKeptCnv <- function(key = NULL, cnv = NULL, sourceId = NULL){
    if(is.null(key)) {
        if(is.null(sourceId)) sourceId <- sourceId()
        key <- getCnvKey(cnv, sourceId)
    }
    !is.null(keptCnvs[[key]]) && keptCnvs[[key]]
}

#----------------------------------------------------------------------
# loade sample/cell source data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
# projectName <- projectNameReactive(sourceId)
project <- normalizeDataReactive(sourceId)
sourceIsInitializing <- TRUE
cells <- reactive({
    project <- project()
    req(project)  
    sourceId <- sourceId()
    projectCellIds <- projectCellIds()
    sourceIsInitializing <<- TRUE
    i <- project$colData[, {
        switch(
            input$cellStatus,
            All    = rep(TRUE, .N),
            Bad    = bad,
            Keep   = getKeep(bad, keep, cell_id, sourceId, FALSE),
            Reject = getRejected(bad, keep, cell_id, sourceId, FALSE)
        ) & switch(
            input$replicating,
            All = TRUE,
            Yes = !bad &  getReplicating(bad, replicating, cell_id, sourceId, FALSE),
            No  = !bad & !getReplicating(bad, replicating, cell_id, sourceId, FALSE)
        ) & (projectCellIds[1] == ALL | cell_id %in% projectCellIds) &
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
    project <- project()
    req(project)    
    sapply(project$colData$cell_id, function(cell_id){
        cell <- project$cells[[cell_id]]
        if(cell$badCell) return(NA)
        shapeModel <- getShapeModel(settings, cell)
        cw <- cell$windows[[shapeModel$key]]
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
# cascade update sample and cell selectors
#----------------------------------------------------------------------
ALL <- "__ALL__"
projectCellIds <- reactiveVal(ALL)
updateCellFilter <- function(){
    project <- project()
    req(project)
    i <- project$manifest[, input$sampleNameFilter == ALL | Sample_Name == input$sampleNameFilter]
    choices <- project$manifest[i, Sample_ID]
    names(choices) <- project$manifest[i, Description]
    projectCellIds(choices)
    freezeReactiveValue(input, "pageNumber")
    updateTextInput(session, "pageNumber", value = 1)
    freezeReactiveValue(input, "cellIdFilter")
    updateSelectInput(session, "cellIdFilter", choices = c(All = ALL, sort(choices)), selected = ALL)    
}
observeEvent(sourceId(), {
    project <- project()
    req(project)
    choices <- project$manifest$Sample_Name
    names(choices) <- choices
    freezeReactiveValue(input, "sampleNameFilter")
    updateSelectInput(session, "sampleNameFilter", choices = c(All = ALL, sort(unique(choices))), selected = ALL)
    updateCellFilter()
})
observeEvent(input$sampleNameFilter, {
    if(sourceIsInitializing) {
        sourceIsInitializing <<- FALSE
        return()
    }
    updateCellFilter()
}, ignoreInit = TRUE)
# observeEvent(cells(), {
#     freezeReactiveValue(input, "pageNumber")
#     updateTextInput(session, "pageNumber", value = 1)
# })

#----------------------------------------------------------------------
# buttons for setting user cell overrides
#----------------------------------------------------------------------
keepRejectButtons <- function(colData, cell){
    prefix <- session$ns("")
    kept <- colData[, getKeep(bad, keep, cell_id, sourceId())]        
    keepRejectOverride <- if(colData$keep){
        if( kept) TRUE else FALSE        
    } else {
        if(!kept) TRUE else FALSE
    }
    if(cell$cellIsReplicating){
        replicating <- getReplicating(cell$badCell, cell$cellIsReplicating, cell$cell_id, sourceId())
        replicatingOverride <- if(replicating) TRUE else FALSE
    }
    tags$div(
        style = "margin-top: 5px;",
        tags$button(
            onclick = paste0("cellToggleOverride('", prefix, "', '", cell$cell_id, "', 'keep', '", keepRejectOverride, "')"),
            if(kept) "Reject" else "Keep"
        ),
        if(cell$cellIsReplicating) tags$button(
            onclick = paste0("cellToggleOverride('", prefix, "', '", cell$cell_id, "', 'replicating', '", replicatingOverride, "')"),
            if(replicating) "Not Replicating" else "Replicating"
        ) else ""
    )
}
observeEvent(input$cellToggleOverride, {
    sourceId <- sourceId()
    x <- input$cellToggleOverride
    overrides[[sourceId]][[x$cell_id]][[x$key]] <- as.logical(x$override)
})

#----------------------------------------------------------------------
# a stack of individual cell plots
#----------------------------------------------------------------------
cellStack <- reactive({   
    cells <- cells()
    if(cells$n == 0) return(integer())
    project <- project()
    sortBy <- settings$get("Page_Options", "Sort_By")   
    ascDesc <- if(settings$get("Page_Options", "Order") == "Ascending") 1 else -1
    colData <- project$colData[cells$i][order(ascDesc * switch(
        sortBy,
        windowPower = project$colData[cells$i, windowPower + bad],
        cnsd = project$colData[cells$i, cnsd],
        fractionS = project$colData[cells$i, fractionS],
        slope = slopes()[cells$i],
        project$colData[cells$i, cell_id]
    ))]
    list(
        colData = colData,
        i = 1:input$cellsPerPage + input$cellsPerPage * (as.integer(input$pageNumber) - 1) 
    )
})
output$genomePlots <- renderUI({     
    cells <- cells()
    req(cells)
    if(cells$n == 0) return("no cells to plot")
    startSpinner(session, message = "loading genome plots")
    cellStack <- cellStack()
    project <- project()
    labelRow <- {
        colData <- cellStack$colData[cellStack$i[1], ]
        cell <- project$cells[[colData$cell_id]]
        createGenomeLabelRow(project, cell)
    }
    cells <- lapply(cellStack$i, function(i){
        if(i > cells$n) return(NULL)
        colData <- cellStack$colData[i, ]
        cell <- project$cells[[colData$cell_id]]
        plotOneCellUI_genome(sourceId(), project, cell, settings, keepRejectButtons, getReplicating, getKeep)
    })
    isolate({ initGenomePlotClicks(initGenomePlotClicks() + 1) })
    stopSpinner(session)
    tagList(labelRow, cells)
})
observeEvent(input$cellsPerPage, {
    session$sendCustomMessage("cellPlotsWrapperUpdate", list(
        prefix = session$ns(""),
        cellsPerPage = input$cellsPerPage
    ))
})

#----------------------------------------------------------------------
# handle clicks into chromosome-level views
#----------------------------------------------------------------------
zoomChrom <- reactiveVal(NULL)
invalidateZoomPlots <- reactiveVal(0)
invalidateZoomedCell <- NULL
initGenomePlotClicks <- reactiveVal(0)
observeEvent(initGenomePlotClicks(), {
    for(divId in c("genomePlotsWrapper", "chromPlotsWrapper")){
        session$sendCustomMessage("cellPlotsWrapperInit", list(
            prefix = session$ns(""),
            divId = divId
        ))
    }
}, ignoreInit = TRUE)
observeEvent(input$cellWindowsPlotClick, {
    handleCellPlotClick(project(), input$cellWindowsPlotClick, zoomChrom, zoomTargetWindow)
})
output$chromPlots <- renderUI({ 
    req(zoomChrom())
    invalidateZoomPlots()
    startSpinner(session, message = "loading chrom plots")
    cells <- cells()
    project <- project()
    cellStack <- cellStack()
    labelRow <- createZoomLabelRow(session, zoomChrom())
    cells <- lapply(cellStack$i, function(i){
        if(i > cells$n) return(NULL)
        colData <- cellStack$colData[i, ]
        cell <- project$cells[[colData$cell_id]]
        plotOneCellUI_chrom(sourceId(), project, cell, settings, zoomChrom(), invalidateZoomedCell, 
                            getReplicating, getKeep, isKeptCnv)
    })
    isolate({ initGenomePlotClicks(initGenomePlotClicks() + 1) })
    invalidateZoomedCell <<- NULL    
    stopSpinner(session)
    tagList(labelRow, cells)
})
prevNextZoomChrom <- function(inc){
    cellStack <- cellStack()
    project <- project()
    zoomChrom <- zoomChrom()
    colData <- cellStack$colData[cellStack$i[1], ]
    cell <- project$cells[[colData$cell_id]]
    chroms <- project$windows[[cell$windowPower + 1]][, unique(chrom)]
    i <- which(chroms == zoomChrom) + inc
    if(i < 1 || i > length(chroms)) return(NULL)
    zoomChrom(chroms[i])
}
observeEvent(input$prevZoomChrom, { prevNextZoomChrom(-1) })
observeEvent(input$nextZoomChrom, { prevNextZoomChrom( 1) })
observeEvent(input$closeZoomPanel, { zoomChrom(NULL) })

#----------------------------------------------------------------------
# handle clicks within chromosome-level views
#----------------------------------------------------------------------
zoomTargetWindow <- reactiveVal(NULL)
observeEvent(zoomTargetWindow(), {
    x <- zoomTargetWindow()
    project <- project()
    sampleName_ <- project$manifest[Sample_ID == x$cell_id, Sample_Name]
    cell <- project$cells[[x$cell_id]]
    shapeModel <- getShapeModel(settings, cell)
    cnvs <- project$cnvs[[shapeModel$key]]
    if(is.null(cnvs) || nrow(cnvs) == 0) return(NULL)
    w <- x$targetWindow
    cnv <- cnvs[
        type == getCnvsType(settings) & 
        sampleName == sampleName_ & 
        cell_id == x$cell_id &
        chrom == w$chrom & 
        start <= w$start & 
        end >= w$end
    ]
    if(nrow(cnv) != 1) return(NULL)
    toggleKeptCnv(cnv, sourceId = NULL)
    invalidateZoomedCell <<- x$cell_id
    invalidateZoomPlots( invalidateZoomPlots() + 1)
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    bo <- bm$outcomes
    settings$replace(bm$settings)
    for(sourceId in names(bo$overrides)) overrides[[sourceId]] <- bo$overrides[[sourceId]]
    updateTextInput(session, 'cellsPerPage', value = bm$input$cellsPerPage)
    updateTextInput(session, 'pageNumber',   value = bm$input$pageNumber)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,  
    outcomes = list(
        overrides = overrides
    ),      
    NULL
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
