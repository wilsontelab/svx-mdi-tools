#----------------------------------------------------------------------
# reactive components to generate interactive plots for optimizing GC bias
# normalization in bin count data and to reject low quality samples/cells
#----------------------------------------------------------------------
# expects a prior appStep such as assignSamples that creates SampleSets
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
normalizeGCServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'normalizeGC' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id
)
sampleSet <- sampleSetGroupTypeServer( # selectors to pick a sample set, group and type
    id = 'data',
    parentId = id
)
stepper <- listStepperButtonsServer( # inputs to pick a single sample or cell from the sample set
    id = 'stepper',
    dataReactive = sampleSet$assignments,
    nameFn = getAssignedSampleName
)
libraryQCFailed <- getStepOutcomesByType('libraryQC')
outcomes <- reactiveValues()

#----------------------------------------------------------------------
# react to sample/cell selectors by collecting ~static source data
#----------------------------------------------------------------------
currentSourceId <- reactiveVal() # only update currentSourceId() and mse() if the source changes
currentSampleId <- reactive({ stepper$getCurrent()$Sample_ID })
observeEvent(stepper$getCurrent(), {
    sample <- stepper$getCurrent()
    req(sample)
    reportProgress(sample$Sample_ID, 'currentSample')
    reportProgress(sample$Source_ID, 'currentSource')
    currentSourceId(sample$Source_ID)
})

# on first encounter of a count matrix, convert it to a raw mdiSummarizedExperiment
# keyed on sourceId only, i.e., mse includes all samples from the source
createMseFromCountMatrix <- function(keyObject, ...){
    reportProgress('create raw mdiSummarizedExperiment', module)
    sourceId <- keyObject
    source <- getSourceFromId(sourceId)
    path <- getSourceFilePath(sourceId, options$contentFileType)   
    new_mdiSummarizedExperiment(path, manifest = source$unique)
}

# on first encounter, or if bin filters change, remove untrustworthy bins and cache the result
# keyed on sourceId and Bin_Filters (but not Fit_Parameters)
createFilteredMse <- function(keyObject, ...){
    reportProgress('create filtered mdiSummarizedExperiment', module)   
    mse <- dataCache$get(CONSTANTS$mseInitialType,
                         keyObject = keyObject$sourceId, createFn = createMseFromCountMatrix, 
                         from = "disk", create = "asNeeded")
    mse <- mse$value # see note below about possible cache purging 
    fs <- keyObject$settings$Bin_Filters
    rowData <- data.table( data.frame( mse$getRowData() ) )
    widths <- getCol(mse, 'width')
    attributes <- names(rowData)       
    rowData[, masked := FALSE] # logical of bins to reject, i.e. mask, i.e. filter against    
    if(fs$Mappability_Source$value %in% attributes){ # bins must be mappable
        mappable <- rowData[[fs$Mappability_Source$value]] >= fs$Minimum_Mappability$value
        rowData[, masked := masked | !mappable] 
    }   
    if('gc' %in% attributes){ # disard bins with extreme outlier GC percentage
        gcOK <- rowData[, gc] >= fs$Minimum_GC$value & rowData[, gc] <= fs$Maximum_GC$value
        rowData[, masked := masked | !gcOK] 
    }         
    if('excluded' %in% attributes){ # discard bins with too few included bases
        rowData[, masked := masked | excluded / widths > fs$Maximum_Excluded$value]
    }     
    rowData[, masked := masked | isLastBin] # discard the last bin in each chromosome, it is untrustworthy
    mse$subset(rows = !rowData[, masked]) 
}

# on first encounter, or if bin filters change, apply pre-optimization calculations to filtered read depth matrix
# keyed on sourceId, Bin_Filters and Fit_Parameters (but as above, includes all samples from a source, even those not in a Sample Set) # nolint
createAssembledMSE <- function(keyObject, ...){
    reportProgress('create pre-optimization mdiSummarizedExperiment', module)
    keyObject$Fit_Parameters <- NULL
    mse <- dataCache$get(CONSTANTS$mseFilteredType,
                         keyObject = keyObject, createFn = createFilteredMse, 
                         from = "disk", create = "asNeeded") # has same settings as filtering step
    mse <- mse$value # see note below about possible cache purging

    # mappability and excluded bases adjustment
    NR_raw_b <- mse$getAssay('raw') # bin read depth
    map_b    <- getCol(mse, 'genmap') # mappability; here the average value of all _non-excluded_ bin positions    
    exc_b    <- getCol(mse, 'excluded') / getCol(mse, 'width')
    gc_b     <- getCol(mse, 'gc') # fraction GC
    NR_map_b <- NR_raw_b / map_b / (1 - exc_b)

    # apply rolling windows if requested
    windowSize <- settings$Fit_Parameters()$Window_Size$value
    if(windowSize > 1) {
        chroms <- getCol(mse, 'chrom')
        NR_map_b[, chrom := chroms]
        gc_b <- data.table(gc_b = gc_b, chrom = chroms)        
        NR_map_b <- NR_map_b[, lapply(.SD, rollsum_permute,  windowSize, na.rm = TRUE), by = 'chrom']
        gc_b <-         gc_b[, lapply(.SD, rollmean_permute, windowSize, na.rm = TRUE), by = 'chrom']
        NR_map_b$chrom <- NULL
        gc_b <- gc_b$gc_b
    }
    
    # suppress 0 double error that report negative numbers (e.g., -0.000000001) mainly due to rollsum_permute on stretches of 0's # nolint
    NR_map_b <- data.table( sapply(seq_len(ncol(NR_map_b)), function(j){ 
        pmax(0L, NR_map_b[[j]], na.rm = TRUE)
    }) )
    setnames(NR_map_b, names(NR_raw_b))
    
    # aggregate NR_map_b; use median instead of mean, more stable to outliers
    NR_map_g <- apply(NR_map_b, 2, median, na.rm = TRUE) 
    mse$setRowData('gc_b_ma', gc_b)

    # initial GC bias correction by negative binomial curve fit to all mappability-adjusted bins
    NR_gc0_b <- data.table( sapply(seq_len(ncol(NR_map_b)), function(j){   
        fit <- new_nbinomCountsGC(NR_map_b[[j]], gc_b, binCN = 1) # CN of 1 means we use the counts as is to begin with
        ER_gc0_b <- predict(fit, gc_b, type = 'adjustedPeak') # use peak for visualization, unless it is unreliable
        pmax(1L, NR_map_b[[j]] * NR_map_g[j] / ER_gc0_b, na.rm = TRUE) # again, suppress negative and zero values of ER_gc0_b # nolint
    }) )
    setnames(NR_gc0_b, names(NR_map_b))

    # calculate a simplified set of Lorenz points for rapid plotting
    # base these on NR_gc0_b; do not worry about updating with NR_gc_b
    counted <- getCol(mse, 'counted')
    lorenzPlotData <- lapply(seq_len(ncol(NR_gc0_b)), function(j){
        y <- NR_gc0_b[[j]]    
        order <- order(y)     
        x <- counted[order]    
        y <- y[order]
        data.table(
            x = cumsum(as.numeric(x)) / sum(x),
            y = cumsum(as.numeric(y)) / sum(y)
        )
    })
    names(lorenzPlotData) <- names(NR_gc0_b)
    lorenzPoints <- seq(1, nrow(NR_gc0_b), length.out = nrow(NR_gc0_b) / 10)
    revLorenzPoints <- rev(lorenzPoints)
    lorenzPlotDataAll <- do.call(rbind, lapply(seq_len(ncol(NR_gc0_b)), function(j){   
        lorenzPlotData[[j]][if(j %% 2 == 0) lorenzPoints else revLorenzPoints]
    }))

    # build the final mdiSummarizedExperiment, keep only the layers we need
    mse$setAssay('NR_map_b', NR_map_b) # used to calculate final NC_gc
    mse$setAssay('NR_gc0_b', NR_gc0_b) # plotted for user selection
    mse$setMetadata('lorenzPlotData', lorenzPlotData)
    mse$setMetadata('lorenzPlotDataAll', lorenzPlotDataAll)
    mse
}

# get the data for a newly selected data source, from cache or building objects as required
mse <- reactive({ # mse = mdiSummarizedExperiment, a version of RangedSummarizedExperiment
    mse <- getGcBiasMse(module, options, currentSourceId(), 
                        create = "asNeeded", createFn = createAssembledMSE)
    mse
})

#----------------------------------------------------------------------
# react to new data load by initializing user-dependent data values, e.g. bin selection
# these are cached separately and also stored in bookmark
#----------------------------------------------------------------------
initializeOutcomes <- function(key, ...){
    if(!is.null(outcomes[[key]])) return( outcomes[[key]]$value ) # return a bookmarked value (not the entire object)
    reportProgress('initializeOutcomes', module)
    mse <- mse()$value
    nSamples  <- nSamples(mse)
    nFeatures <- nFeatures(mse)
    sampleIds <- colnames(mse$rse())
    CN_b_rle <- rle(rep(CONSTANTS$noCnCall, nFeatures)) # therefore, all samples require user to select at least some modal CN bins (or reject the sample) # nolint
    x <- list(
        batchReject_s = as.list(rep(CONSTANTS$batched, nSamples)), # whether a sample has been rejected (by default they are fully included) # nolint
        MCN_s = as.list(rep(CONSTANTS$noCnCall, nSamples)), # sample-specific modal CN; default triggers an auto-fill from the project default MCN on first load # nolint
        CN_b = rep(list(CN_b_rle), nSamples), # the bin spans assigned by user as having a specific CN state # nolint
        fit = as.list(rep(NA, nSamples)), # the revised negative binomial model calculated by nbinomCountsGC based on CN_b # nolint
        hmm = as.list(rep(NA, nSamples))  # the optimized HMM fit to CN calculated from fit, with user adjustments
    ) 
    for(name in c('batchReject_s', 'MCN_s', 'CN_b', 'fit', 'hmm')) names(x[[name]]) <- sampleIds
    x
}
observeEvent(mse(), {
    req(mse())
    startSpinner(session, paste(module, 'get outcomes list'))
    # NB: key, not cacheKey, we match outcomes keys to mse keys, i.e., static and outcomes objects are keyed the same
    # thus, a sample will behave exactly the same in two different sample sets; editing one place edits it every place
    outcomes[[mse()$key]] <- dataCache$get(CONSTANTS$gcBiasOutcomesType, keyObject = mse()$keyObject,
                                           from = "ram", permanent = FALSE,
                                           create = "asNeeded", createFn = initializeOutcomes)
    stopSpinner(session) 
})

#----------------------------------------------------------------------
# sample outcome getters and setters (setSampleOutcome is the ultimate point of this UI)
#----------------------------------------------------------------------
sourceInfo <- reactive({
    d <- sampleSet$assignments()     
    req(d)
    req(nrow(d) > 0)
    i <- stepper$current()    
    req(i)
    maxI <- nrow(d)
    if(i > maxI) i <- 1
    x <- d[i, ]
    list(
        sourceId = x$Source_ID,
        sampleId = x$Sample_ID,
        failed   = libraryQCFailed()[[x$Source_ID]][x$Sample_ID]
    ) 
})
getSampleOutcome <- function(type, allSamples = FALSE){
    if(allSamples) return( outcomes[[mse()$key]]$value[[type]] )
    else {
        si <- sourceInfo()
        req(si)
        outcomes[[mse()$key]]$value[[type]][[si$sampleId]]        
    }
}
setSampleOutcome <- function(type, value){
    si <- sourceInfo()    
    req(si)
    outcomes[[mse()$key]]$value[[type]][[si$sampleId]] <- value
    dataCache$set(outcomes[[mse()$key]])
    updateSampleTableBuffer(si$sampleId, type, value)
    if(type == 'CN_b') calculateRevisedFit(si$sampleId, value)
    value # echo back the newly commited value
}

#----------------------------------------------------------------------
# establish keys for invalidating plot images
#----------------------------------------------------------------------
preOutcomesSampleKey <- reactive({ # depend on the our input summarized experiment and selected sample
    list(                          # but not any user interactive outcomes
        mseKey = mse()$key,
        source = sourceInfo(),
        failed = libraryQCFailed()
    )
})
sampleOutcomesKey <- reactive({ # additionally depend on the user's interactive outcomes
    x <- preOutcomesSampleKey()
    x$outcomes <- lapply(c('MCN_s', 'CN_b', 'fit'), function(type) {
        outcomes[[x$mseKey]]$value[[type]][[x$source$sampleId]]
    })
    x
})

#----------------------------------------------------------------------
# getting and setting of MCN_s and batchReject_s
#----------------------------------------------------------------------
# update sample-level inputs when user changes samples
observe({
    req(mse())
    req(stepper$current())
    req(sampleSet$assignments())
    isolate({
        mcn <- getSampleOutcome('MCN_s')
        if(mcn == CONSTANTS$noCnCall) mcn <- setSampleOutcome('MCN_s', settings$Fit_Parameters()$Default_Modal_CN$value)
        updateNumericInput(session, 'sampleModalCN', value = mcn)
        updateRadioButtons(session, 'batchRejectSample',  selected = getSampleOutcome('batchReject_s'))    
    })
})
# react to sample-level outcome changes by user (modal CN, sample rejection)
observeEvent(input$sampleModalCN, {
    req(input$sampleModalCN)
    MCN_s_prev <- getSampleOutcome('MCN_s')
    MCN_s_new  <- input$sampleModalCN # value is integer
    req(MCN_s_prev != MCN_s_new)
    CN_b <- getSampleOutcome('CN_b')
    assigned <- CN_b$values != CONSTANTS$noCnCall
    CN_b$values[assigned] <- CN_b$values[assigned] + (MCN_s_new - MCN_s_prev)
    setSampleOutcome('CN_b', CN_b)
    setSampleOutcome('MCN_s', MCN_s_new)
})
observeEvent(input$batchRejectSample, {
    req(!is.null(input$batchRejectSample))
    setSampleOutcome('batchReject_s', as.integer(input$batchRejectSample)) 
})

#----------------------------------------------------------------------
# react to NR_gc0_b genomic plot selection by setting CN_b
#----------------------------------------------------------------------
observeEvent(gc0ReadDepthPlot$selected(), {
    selected <- gc0ReadDepthPlot$selected()
    req(selected)
    req(length(selected$pointNumber) > 2)
    reportProgress('observe gc0ReadDepthPlot$selected()')
    
    # update CN_b
    mse <- mse()$value
    nFeatures <- nFeatures(mse)
    MCN_b <- 1:nFeatures %in% (selected$pointNumber + 1) # plot_ly is 0-indexed 
    MCN_s <- getSampleOutcome('MCN_s')    
    CN_b <- inverse.rle(getSampleOutcome('CN_b'))
    CN_b[] <- CONSTANTS$noCnCall # MCN selection necessarily clears all prior state specifications
    CN_b[MCN_b] <- MCN_s
    setSampleOutcome('CN_b', rle(CN_b))
})

#----------------------------------------------------------------------
# react to CN_gc genomic plot click by forcing CN_b for a clicked bin
#----------------------------------------------------------------------
observeEvent(gcReadDepthPlot$clicked(), {
    clicked <- gcReadDepthPlot$clicked()
    req(clicked)
    req(clicked$curveNumber == 0) # must have clicked a data point (not the HMM trace)
    binN <- clicked$x       
    CN <- round(clicked$y, 0)       
    CN_b <- inverse.rle(getSampleOutcome('CN_b'))
    CN_b[binN] <- CN    
    setSampleOutcome('CN_b', rle(CN_b))
})

#----------------------------------------------------------------------
# cascade to calculate revised negative binomial fit whenever CN_b changes; stored as an outcome
#----------------------------------------------------------------------
calculateRevisedFit <- function(sampleId, CN_b){
    reportProgress('calculateRevisedFit')
    CN_b <- inverse.rle(CN_b)
    req(sum(CN_b > 0) > 2)
    mse <- mse()$value
    req(mse)
    gc_b     <- getCol(mse, 'gc_b_ma')
    NR_map_b <- getCol(mse, sampleId, 'NR_map_b')

    # fit once to get model based on selected bins only
    nb  <- new_nbinomCountsGC(NR_map_b, gc_b, binCN = CN_b)
    hmm <- viterbi(nb, NR_map_b, gc_b, forceCNs = CN_b, asRle = FALSE)

    # then refit to all bins based on HMM initial estimate, while respecting user-forced bins
    cn_estimate <- ifelse(hmm$cn == hmm$maxCN, NA, hmm$cn) # bins at maxCN are unreliable as they might be >maxCN
    nb  <- new_nbinomCountsGC(NR_map_b, gc_b, binCN = cn_estimate) # this line used the HMM-updated CN for all bins
    hmm <- viterbi(nb, NR_map_b, gc_b, forceCNs = CN_b) # NB: this line still used CN_b as forceCNs (but nb is updated) 
    setSampleOutcome('fit', nb)
    setSampleOutcome('hmm', hmm$cn)
}

#----------------------------------------------------------------------
# cascade to calculate revised CN_gc for plotting; use distribution _peak_ for better visual comparisons
# not stored in outcomes object, but is cached per session for faster sample switching
#----------------------------------------------------------------------
CN_gc <- reactive({
    fit <- getSampleOutcome('fit')
    req(fit)
    req(is.list(fit))
    mse <- mse()$value
    gc_b <- getCol(mse, 'gc_b_ma')
    req(gc_b)
    RPA <- predict(fit, gc_b, type = 'adjustedPeak')
    sourceInfo <- sourceInfo()
    req(sourceInfo)
    req(!sourceInfo$failed)
    getCol(mse, sourceInfo$sampleId, 'NR_map_b') / RPA    
})
MCN_b <- reactive({
    CN_b <- inverse.rle(getSampleOutcome('CN_b'))
    CN_b == getSampleOutcome('MCN_s')    
})

#----------------------------------------------------------------------
# show failed sample feedback
#----------------------------------------------------------------------
failedQCText <- reactive({
    sourceInfo <- sourceInfo()
    req(sourceInfo)
    req(sourceInfo$failed)
    "Failed Library QC"    
}) 
output$failedQC1 <- renderText({ failedQCText() })
output$failedQC2 <- renderText({ failedQCText() })

#----------------------------------------------------------------------
# plot raw and normalized result by genomic bin; enable selection of modal CN bins
#----------------------------------------------------------------------
range_noChrY <- function(d, axis){
    d <- d[getCol(mse()$value, 'chrom') != 'chrY']
    range_both(d, axis, 3)
}
gc0ReadDepthPlotData <- reactive({
    req(mse())
    mse <- mse()$value
    sourceInfo <- sourceInfo()
    req(sourceInfo)
    req(!sourceInfo$failed)
    data.table(
        x = getCol(mse, 'binN'),    
        y = log10( getCol(mse, sourceInfo$sampleId, 'NR_gc0_b') )
    )    
})
gc0ReadDepthPlot <- interactiveScatterplotServer(
    'gc0ReadDepthPlot',
    gc0ReadDepthPlotData,
    accelerate = TRUE,
    xtitle = "Genomic Bin",
    xrange = range_pad,
    xzeroline = FALSE,
    ytitle = 'log10(# Reads)',
    yrange = range_noChrY, 
    selectable = 'h',
    cacheReactive = preOutcomesSampleKey
)
#----------------------------------------------------------------------
getChromLines <- function(d){
    L <- d$x[ getCol(mse()$value, 'isLastBin') ]
    structure(c(0, L) + 0.5, color = rep('grey', length(L) + 1))    
}
getCNLines <- function(...){
    structure(input$sampleModalCN + -2:2,
              color = c('grey', 'grey', 'black', 'grey', 'grey'))   
}
gcReadDepthPlotAllPoints <- reactive({
    req(mse())
    mse <- mse()$value
    CN_gc <- CN_gc()
    req(CN_gc)
    data.table(
        x = getCol(mse, 'binN'),
        y = CN_gc
    )    
})
gcReadDepthOverplotData <- reactive({
    d <- copy(gcReadDepthPlotAllPoints())
    CN_hmm <- getSampleOutcome('hmm')
    req(CN_hmm)      
    d[, ':='(x = x, y = inverse.rle(CN_hmm))]
})
gcReadDepthPlot <- interactiveScatterplotServer(
    'gcReadDepthPlot',
    gcReadDepthPlotAllPoints,
    accelerate = TRUE,
    xtitle = "Genomic Bin",
    xrange = range_pad,
    xzeroline = FALSE,
    ytitle = 'Copy Number',
    yrange = c(0, 5), # TODO: expose maxCN as plotting option
    clickable = TRUE,
    grid = list(x = TRUE, y = FALSE),
    hLines = getCNLines,
    vLines = getChromLines,
    overplot = gcReadDepthOverplotData,
    cacheReactive = sampleOutcomesKey 
)                             

#----------------------------------------------------------------------
# QC plot: distribution of CN_gc from selected/unselected bins
#----------------------------------------------------------------------
yrange_CN_density <- function(...){
    selected   <- gcReadDepthDensityPlotData_MCN()
    unselected <- gcReadDepthDensityPlotData_CNV()
    maxSelected <- max(selected$y, na.rm = TRUE)
    maxUnselected <- if(is.data.frame(unselected)){
        if(nrow(unselected) == 0) 0 else max(unselected$y, na.rm = TRUE)
    } else {
        sapply(unselected, function(x) max(x$y, na.rm = TRUE))
    }
    c(0, max(maxSelected, maxUnselected, na.rm = TRUE)) * 1.1
}
gcReadDepthDensityPlotData_ <- function(negate){
    nullData <- data.table(x = integer(), y = double())
    CN_gc <- CN_gc()
    req(CN_gc)
    CN_hmm <- getSampleOutcome('hmm')
    req(CN_hmm)      
    CN_hmm <- inverse.rle(CN_hmm)
    mcn <- input$sampleModalCN    
    getDT <- function(bins){
        d <- density(CN_gc[bins])
        data.table(x = d$x, y = d$y)            
    }
    if(negate){     
        cns <- sort(unique(CN_hmm))
        cns <- cns[cns %notin% c(mcn, CONSTANTS$noCnCall, 0)]
        if(length(cns) <= 0) return( nullData )
        x <- lapply(cns, function(cn) getDT(CN_hmm == cn))
        names(x) <- cns
        x
    } else {
        getDT(CN_hmm == mcn)
    }  
}
gcReadDepthDensityPlotData_MCN <- reactive({
    gcReadDepthDensityPlotData_(FALSE)
})
gcReadDepthDensityPlotData_CNV <- reactive({
    gcReadDepthDensityPlotData_(TRUE)
})
interactiveScatterplotServer(
    'gcReadDepthDensityPlot',
    gcReadDepthDensityPlotData_MCN, 
    xtitle = "Copy Number",
    ytitle = "Density",
    xrange = c(0, 5),
    yrange = yrange_CN_density,
    grid = list(x = FALSE, y = TRUE),
    vLines = getCNLines,
    mode = 'lines',
    overplot = gcReadDepthDensityPlotData_CNV,
    cacheReactive = sampleOutcomesKey 
)

#----------------------------------------------------------------------
# QC plot: NR_map_b vs. GC, to visualize quality of GC bias fit
#----------------------------------------------------------------------
gcBiasPlotData <- reactive({ # all data points on the GC bias plot, i.e. all bins
    req(mse())
    sourceInfo <- sourceInfo()
    req(sourceInfo)
    req(!sourceInfo$failed)  
    getXY(mse()$value,
          xcol = 'gc',
          ycol = sourceInfo$sampleId,
          assay = 'NR_map_b') # use mappability corrected counts
})
gcBiasPlotFit <- reactive({ 
    fit <- getSampleOutcome('fit')
    if(is.na(fit)) return(NULL)
    data.table(
        x = fit$model$fractionGC,
        y = predict(fit, type = 'peak') * input$sampleModalCN # again, use _peak_ for more intuitive visualization
    )                                                         # here, do _not_ use adjustedPeak
})
gcBiasPlotHighlight <- reactive({ # only the bins selected by the user
    d <- gcBiasPlotData()
    req(d)
    MCN_b <- MCN_b()
    if(sum( MCN_b) == 0 || # if no bins selected
       sum(!MCN_b) == 0){  # or all bins selected
        return( d[FALSE] ) # then no overplot is needed
    }
    d[MCN_b]
})
readDepthByGCPlot <- interactiveScatterplotServer(
    'readDepthByGCPlot',
    gcBiasPlotData,
    accelerate = TRUE,
    xtitle = "Fraction GC",
    ytitle = "# Reads",
    xrange = c(settings$Bin_Filters()$Minimum_GC,
               settings$Bin_Filters()$Maximum_GC),
    yrange = range_pos,
    fitMethod = gcBiasPlotFit,
    overplot = gcBiasPlotHighlight,
    cacheReactive = sampleOutcomesKey 
)

#----------------------------------------------------------------------
# QC plot: sample read depths, to visualize relative sample coverage
#----------------------------------------------------------------------
sampleReadDepths <- reactive({
    req(mse())
    colData <- data.table( data.frame( mse()$value$getColData() ) )
    batchReject_s <- getSampleOutcome('batchReject_s', allSamples = TRUE)
    alpha <- sapply(getSampleOutcome('fit', allSamples = TRUE), function(fit){
        if(is.list(fit)) 1 / fit$theta else -0.1
    })
    data.table(
        Sample_ID = colData$Sample_ID,
        x = colData$sum, # i.e. sample read depth
        y = alpha,
        pointSize = ifelse(colData$Sample_ID %in% sampleSet$assignments()$Sample_ID, 7, 4),
        symbol = ifelse(batchReject_s == CONSTANTS$rejected, 'x', 'circle')
    )
})
sampleReadDepthsPlot <- interactiveScatterplotServer(
    'sampleReadDepthsPlot',
    sampleReadDepths,
    clickable = TRUE, # allow moving to a sample in the current group/type by point click
    xtitle = "Total Read Depth",
    ytitle = "Overdispersion Alpha",
    yrange = function(d, ...) c(-0.15, max(max(0.25, d$y, na.rm = TRUE) + 0.05)),
    overplotPointSize = 7,
    overplot = function(d) d[Sample_ID == currentSampleId()],
    cacheReactive = sampleOutcomesKey 
)
observeEvent(sampleReadDepthsPlot$clicked(), {
    d <- sampleReadDepthsPlot$clicked()
    req(d)
    if(d$curveNumber != 0) return(NULL) # user clicked the point for the current sample
    Sample_ID <- sampleReadDepths()[d$pointNumber + 1, Sample_ID]
    allowed <- sampleSet$assignments()$Sample_ID
    if(Sample_ID %notin% allowed) return(NULL) # user clicked a sample not in the current Group/Type
    stepper$setCurrent(which(allowed == Sample_ID))
})

#----------------------------------------------------------------------
# QC plot: Lorenz plot, to visualize relative sample quality
#----------------------------------------------------------------------
lorenzPlotData <- reactive({
    req(mse())
    mse()$value$getMetadata()$lorenzPlotDataAll
})
interactiveScatterplotServer(
    'lorenzPlot',
    lorenzPlotData,
    xtitle = "Fraction of Genome",
    ytitle = "Fraction of Reads",
    mode = 'lines',
    unityLine = TRUE,
    color = "lightgrey",
    lineWidth = 2,
    overplotLineWidth = 3,
    overplot = function(d){
        req(mse())
        mse()$value$getMetadata()$lorenzPlotData[[ currentSampleId() ]]
    },
    cacheReactive = preOutcomesSampleKey 
)

#----------------------------------------------------------------------
# summary table of all samples, react to table click
#----------------------------------------------------------------------
samplesTable <- bufferedTableServer(
    id = 'samplesTable',
    parentId = id,
    parentInput = input,
    selection = 'single',
    selectionFn = stepper$setCurrent,
    options = list(paging = FALSE),
    tableData = reactive({
        ssa <- sampleSet$assignments() # thus, invalidates only when sample set+group+type changes
        req( outcomes[[mse()$key]] )
        isolate({
            batchReject_s <- getSampleOutcome('batchReject_s', allSamples = TRUE)
            status <- getSampleFlagStatus(unlist(batchReject_s))
            failed <- libraryQCFailed()        
            getSampleFailed <- Vectorize(function(Source_ID, Sample_ID) failed[[Source_ID]][Sample_ID])
            failedQC <- 'Failed QC'
            CN_b <- getSampleOutcome('CN_b', allSamples = TRUE)
            getSampleFitted <- Vectorize(function(Sample_ID) any(CN_b[[Sample_ID]]$values != CONSTANTS$noCnCall))
            fits <- getSampleOutcome('fit', allSamples = TRUE)
            MCN_s <- getSampleOutcome('MCN_s', allSamples = TRUE)
            getSampleAlpha <- Vectorize(function(Sample_ID){
                if(is.list(fits[[Sample_ID]])) round(1 / fits[[Sample_ID]]$theta, 3) else '-'
            })
            dt <- ssa[, .SD, .SDcols = c('Source_ID', 'Project', 'Sample_ID')]
            dt[, ':='(
                sampleName = getSampleNames(sampleIds = Sample_ID),
                modalCN = unlist(MCN_s[Sample_ID]),
                qcStatus = ifelse(getSampleFailed(Source_ID, Sample_ID), failedQC, status[Sample_ID])
            )]
            dt[, ':='(
                modalCN = ifelse(modalCN < 0, '', modalCN),
                fitStatus = ifelse(qcStatus == failedQC, 
                                   failedQC, 
                                   ifelse(getSampleFitted(Sample_ID), 'Fitted', 'Pending')),
                alpha = getSampleAlpha(Sample_ID)
            )]
            dt$Source_ID <- NULL
            dt             
        })
    })
)
observeEvent(stepper$current(), {
    samplesTable$selectRows(stepper$current())
})
updateSampleTableBuffer <- function(sampleId, type, value){
    updateCell <- function(col, value) samplesTable$updateCell(sampleId, col, value, 'Sample_ID')
    if(type == 'MCN_s'){
        updateCell('modalCN', value)
    } else if(type == 'batchReject_s'){
        updateCell('qcStatus', getSampleFlagStatus(value))
    } else if(type == 'fit' && !is.na(value)){
        updateCell('fitStatus', 'Fitted')
        updateCell('alpha', round(1 / value$theta, 3))
    }
}

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    if(!is.null(bm$bookmarkedOutcomes)) bm$outcomes <- bm$bookmarkedOutcomes # legacy name change
    if(!is.null(bm$outcomes)) outcomes <<- listToReactiveValues(bm$outcomes)
    playbackStepBookmark(id, module, session, bookmark, handlers = list(
        'stepper-current' = stepper$overrideDefault, # handling special cascading events in replay
        'data-group' = sampleSet$overrideDefault,
        'data-type'  = sampleSet$overrideDefault
    ))
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = reactive({ reactiveValuesToList(outcomes) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
