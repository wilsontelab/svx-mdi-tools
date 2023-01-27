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
isUserOverride <- Vectorize(function(cell_id, key, sourceId) {
    !is.null(overrides[[sourceId]][[cell_id]][[key]]) && overrides[[sourceId]][[cell_id]][[key]]
})
getKeep <- function(bad, keep, cell_id, sourceId = NULL) {
    if(is.null(sourceId)) sourceId <- sourceId()
    overridden <- isUserOverride(cell_id, "keep", sourceId)
    (!bad &  keep & !overridden) | 
    (!bad & !keep &  overridden)
}
getRejected <- function(bad, keep, cell_id, sourceId = NULL){
    if(is.null(sourceId)) sourceId <- sourceId()
    overridden <- isUserOverride(cell_id, "keep", sourceId)
    (!bad & !keep & !overridden) | 
    (!bad &  keep &  overridden)
}
getReplicating <- function(bad, replicating, cell_id, sourceId = NULL){
    if(is.null(sourceId)) sourceId <- sourceId()
    overridden <- isUserOverride(cell_id, "replicating", sourceId)
    !bad & replicating & !overridden
}

#----------------------------------------------------------------------
# source data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
# projectName <- projectNameReactive(sourceId)
project <- normalizeDataReactive(sourceId)
cells <- reactive({
    project <- project()
    req(project)  
    sourceId <- sourceId()
    projectCellIds <- projectCellIds()
    i <- project$colData[, {
        switch(
            input$cellStatus,
            All    = rep(TRUE, .N),
            Bad    = bad,
            Keep   = getKeep(bad, keep, cell_id, sourceId),
            Reject = getRejected(bad, keep, cell_id, sourceId)
        ) & switch(
            input$replicating,
            All = TRUE,
            Yes = !bad &  getReplicating(bad, replicating, cell_id, sourceId),
            No  = !bad & !getReplicating(bad, replicating, cell_id, sourceId)
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
        shapeModel <- settings$get("Page_Options", "Shape_Model")
        cw <- cell$windows[[shapeModel]]
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
    updateCellFilter()
})
observeEvent(cells(), {
    freezeReactiveValue(input, "pageNumber")
    updateTextInput(session, "pageNumber", value = 1)
})

#----------------------------------------------------------------------
# a stack of individual cell plots
#----------------------------------------------------------------------
output$cellPlots <- renderUI({ 
    sourceId <- sourceId()     
    cells <- cells()
    req(cells)
    if(cells$n == 0) return("no cells to plot")
    startSpinner(session, message = "loading cell plots")
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
    cellIndices <- 1:input$cellsPerPage + input$cellsPerPage * (as.integer(input$pageNumber) - 1)
    x <- lapply(cellIndices, function(i){
        if(i > cells$n) return(NULL)
        colData <- colData[i, ]
        cell <- project$cells[[colData$cell_id]]
        plotOneCellUI(project, cell, settings, keepRejectButtons, getReplicating)
    })
    stopSpinner(session)
    x
})

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
