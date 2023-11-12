#----------------------------------------------------------------------
# server components for the normalizeGC appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
normalizeGCServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'normalizeGC'
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
gcPlotModels <- reactiveValues()

#----------------------------------------------------------------------
# sample sources
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
) 
sampleNames <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    source <- getSourceFromId(sourceId)
    req(source)
    source$manifest$Sample_ID
})
coverageTable <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    sampleNames <- sampleNames()
    req(sampleNames)
    startSpinner(session, message = "loading coverage")
    coverageFile <- getSourceFilePath(sourceId, "coverageFile")
    readRDS(coverageFile)[, .SD, .SDcols = c("gc", sampleNames)]    
})
cnvCoverageTable <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading CNV coverage")
    cnvCoverageFile <- getSourceFilePath(sourceId, "cnvCoverageFile")
    x <- readRDS(cnvCoverageFile)
    setkey(x, svId)
    x
})
cnvs <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading CNVs")
    svsFile <- getSourceFilePath(sourceId, "structuralVariants")
    x <- readRDS(svsFile)[JXN_TYPE %in% c("D", "L"), .(
        svId = SV_ID,
        JXN_TYPE,
        SV_SIZE,
        N_SAMPLES
    )]
    setkey(x, svId)
    x
})

#----------------------------------------------------------------------
# interactive GC bias plot, selection cascades to solving negative binomial
#----------------------------------------------------------------------
gcPlotObservers <- list()
gcPlotData <- function(sampleName){
    x <- coverageTable()   
    startSpinner(session, message = paste("plotting", sampleName))
    x <- x[, .SD, .SDcols = c("gc", sampleName)]   
    setnames(x, c("x", "y"))
    x[sample.int(.N, min(.N, 5000))]
}
gcOverplotData <- function(sourceId, sampleName){
    fit <- gcPlotModels[[sourceId]][[sampleName]]  
    if(!isTruthy(fit)) return(NULL)
    startSpinner(session, message = paste("overplotting", sampleName))
    gc <- fit$model$fractionGC
    rpa <- predict(fit, gc, type = 'mu')
    data.table(
        x = c(gc, NA, gc, NA, gc),
        y = c(rpa * 1, NA, rpa * 2, NA, rpa * 3)
    )
}
gcPlotServer <- function(id, sourceId, sampleName){
    plot <- interactiveScatterplotServer(
        id,
        plotData = reactive({ 
            x <- gcPlotData(sampleName) 
            stopSpinner(session)
            x
        }),
        accelerate = TRUE,
        color = CONSTANTS$plotlyColors$blue,
        overplot = reactive({
            x <- gcOverplotData(sourceId, sampleName)
            stopSpinner(session)
            x
        }),
        overplotMode = "lines",
        overplotColor = CONSTANTS$plotlyColors$red,
        xtitle = "Fraction GC",
        xrange = c(0.2, 0.75),
        ytitle = "Read Depth",
        yrange = function(...) range_pos(..., foldIQR = 3),
        selectable = "lasso"
    )    
    if(!is.null(gcPlotObservers[[sampleName]])) gcPlotObservers[[sampleName]]$destroy()
    gcPlotObservers[[sampleName]] <<- observeEvent(plot$selected(), {
        selected <- plot$selected()
        req(selected, nrow(selected) > 10)
        selected <- as.data.table(selected)
        startSpinner(session, message = paste("fitting", sampleName, "GC"))
        if(is.null(gcPlotModels[[sourceId]])) gcPlotModels[[sourceId]] <- list()
        gcPlotModels[[sourceId]][[sampleName]] <- new_nbinomCountsGC(binCounts = selected$y, fractionGC = selected$x, method = 'cubic')
        stopSpinner(session)
    })
}
gcPlotUI <- function(id){
    column(
        width = 6,
        interactiveScatterplotUI(session$ns(id), height = '400px')
    )    
}

#----------------------------------------------------------------------
# interactive CNV plot, correlating CNV sizes to normalized CN
#----------------------------------------------------------------------
cnvPlotData <- function(sourceId, sampleName){
    startSpinner(session, message = paste("normalizing", sampleName, "CNVs"))
    d <- getCnvNormalizedCN(sourceId, sampleName)
    d <- d[, .(
        x = log10(SV_SIZE),
        y = cn,
        JXN_TYPE
    )]
    list(
        Del = d[JXN_TYPE == "L"],        
        Dup = d[JXN_TYPE == "D"]
    )
}
cnvPlotServer <- function(id, sourceId, sampleName){
    plot <- interactiveScatterplotServer(
        id,
        plotData = reactive({ 
            x <- cnvPlotData(sourceId, sampleName) 
            stopSpinner(session)
            x
        }),
        accelerate = TRUE,
        # symbol = "color", #CONSTANTS$plotlyColors$blue,
        pointSize = 5,
        xtitle = "log10(CNV Size)",
        xrange = c(2, 7),
        ytitle = "Copy Number",
        yrange = c(0, 5)
    ) 
}
cnvPlotUI <- function(id){
    column(
        width = 6,
        interactiveScatterplotUI(session$ns(id), height = '400px')
    )    
}

#----------------------------------------------------------------------
# one normalization box per sample, with GC and CNV plots
#----------------------------------------------------------------------
output$samples <- renderUI({
    sourceId <- sourceId()
    sampleNames <- sampleNames()
    req(sourceId, sampleNames)
    lapply(sampleNames, function(sampleName){
        gcBiasId  <- paste("gcBias",  sampleName, sep = "_")
        cnvPlotId <- paste("cnvPlot", sampleName, sep = "_")
        gcPlotServer( gcBiasId, sourceId, sampleName)
        cnvPlotServer(cnvPlotId, sourceId, sampleName)
        fluidRow(box(
            width = 12,
            title = sampleName,
            solidHeader = TRUE,
            status = "primary",
            gcPlotUI(gcBiasId),
            cnvPlotUI(cnvPlotId)
        ))
    })
})

#----------------------------------------------------------------------
# normalization functions
#----------------------------------------------------------------------
getBinNormalizedCN <- function(sourceId, sampleName, gc, coverage){
    fit <- gcPlotModels[[sourceId]][[sampleName]]  
    if(!isTruthy(fit)) return(NULL)
    allowedGC <- range(fit$model$fractionGC)
    rpa <- predict(fit, pmax(allowedGC[1], pmin(allowedGC[2], gc)), type = 'mu')
    coverage / rpa
}
getCnvNormalizedCN <- function(sourceId, sampleName, svIds = NULL){
    fit <- gcPlotModels[[sourceId]][[sampleName]]  
    cnvCoverage <- cnvCoverageTable()[sample == sampleName]
    cnvs <- cnvs()
    req(fit, cnvCoverage, cnvs) 
    if(is.null(svIds)) {
        uniqueSvIds <- cnvs[N_SAMPLES == 1, svId]
        svIds <- cnvCoverage[svId %in% uniqueSvIds, svId] # [1:500]
    }
    allowedGC <- range(fit$model$fractionGC)
    merge(
        cnvCoverage[svIds, {
            gc       <- as.numeric(strsplit(gc, ",")[[1]])
            coverage <- as.numeric(strsplit(coverage, ",")[[1]])
            overlap  <- as.integer(strsplit(overlap, ",")[[1]])
            rpa <- predict(fit, pmax(allowedGC[1], pmin(allowedGC[2], gc)), type = 'mu')
            cn <- coverage / rpa
            .(cn = weighted.mean(cn, overlap))
        }, by = .(svId)],
        cnvs,
        by = "svId",
        all.x = TRUE
    )
}

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    for(sourceId in names(bm$outcomes$gcPlotModels)){
        gcPlotModels[[sourceId]] <- bm$outcomes$gcPlotModels[[sourceId]]
    }
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
    outcomes = list(
        gcPlotModels = gcPlotModels
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    getBinNormalizedCN = getBinNormalizedCN,
    getCnvNormalizedCN = getCnvNormalizedCN,
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
