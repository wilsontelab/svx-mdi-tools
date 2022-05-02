#----------------------------------------------------------------------
# reactive components to generate stacked genomic plots and heatmaps
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
compositePlotsServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'compositePlots' # for reportProgress tracing
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
gcBiasSettings <- getStepSettingsByType("normalizeGC")

#----------------------------------------------------------------------
# react to sample/source selectors by collecting data from GC bias and batch effects normalization
#----------------------------------------------------------------------
sourceInfo <- reactive({
    assignments <- sampleSet$assignments()
    req(assignments)
    list(
        sourceId  = assignments[1, Source_ID],
        sampleIds = assignments[,  Sample_ID],
        key = paste(assignments[,  Sample_ID], collapse = "~")
    )
})
mse <- reactive({ # mse = mdiSummarizedExperiment, a version of RangedSummarizedExperiment
    sourceInfo <- sourceInfo()
    req(sourceInfo)
    getGcBiasMse(module, options, sourceInfo$sourceId,
                 sampleIds = sampleSet$assignments()$Sample_ID)
})
gcBias <- reactive({
    getGcBiasOutcomes(mse)
})
batchEffects <- reactive({
    getBatchEffects(module, options, mse(), gcBias(), sampleSet$allAssignments()$Sample_ID)
})

#----------------------------------------------------------------------
# cascade to calculate copy number per bin and hierarchically cluster based on CN
#----------------------------------------------------------------------
calculateCopyNumber <- function(keyObject, batchEffects, ...){
    reportProgress('calculateCopyNumber', module)
    mse <- mse()$value    
    CN_be <- batchEffects$value$readPerAllele_peak # variable named for the eventual contents, not RPA  
    for(sampleId in names(CN_be)){                 # use peak for more natural visualization
       CN_be[[sampleId]] <-  getCol(mse, sampleId, 'NR_map_b') / CN_be[[sampleId]]
    }
    CN_be[, ':='(
        binN = 1:.N,
        isLastBin = getCol(mse, 'isLastBin')
    )]
    CN_be
}
copyNumber <- reactive({
    batchEffects <- batchEffects()
    req(batchEffects)
    startSpinner(session, paste(module, 'getting copy number'))
    cn <- dataCache$get(CONSTANTS$compositeCopyNumberType, keyObject = batchEffects$keyObject,
                        from = "ram", permanent = TRUE,
                        create = "asNeeded", createFn = calculateCopyNumber,
                        batchEffects = batchEffects)
    #stopSpinner(session)
    req(cn)
    cn$value    
})

#----------------------------------------------------------------------
# sample listing and sorting
#----------------------------------------------------------------------
goodSourceSamples <- reactive({
    batchEffects <- batchEffects()
    req(batchEffects)    
    intersect(batchEffects$value$sampleIds,
              sampleSet$assignments()$Sample_ID)
})
MCN_s <- reactive({
    sampleIds <- goodSourceSamples()
    unlist(gcBias()$MCN_s[sampleIds])
})
clusteredSamples <- reactive({
    copyNumber <- copyNumber()
    req(copyNumber)
    sampleIds <- goodSourceSamples()
    cns <- as.matrix(copyNumber[, .SD, .SDcols = sampleIds])
    MCN_s <- MCN_s()
    cncs <- sapply(seq_len(ncol(cns)), function(j) cns[, j] - MCN_s[j])
    dist <- dist(t( abs(cncs) )) # TODO: enable settings to determine clustering
    hc <- hclust(dist)
    sampleIds[hc$order]
})

#----------------------------------------------------------------------
# cascade to solve HMM for each sample
#----------------------------------------------------------------------
createHMMs <- function(keyObject, gcBias, batchEffects, ...){
    reportProgress('createHMMs', module)
    mse <- mse()$value
    gc_b <- getCol(mse, 'gc_b_ma')
    percentiles <- batchEffects$value$lowertail$median
    sampleIds <- names(batchEffects$value$readPerAllele_peak)
    hmm <- lapply(sampleIds, function(sampleId){
        viterbi(
            gcBias$fit[[sampleId]],
            binCounts = getCol(mse, sampleId, 'NR_map_b'),
            fractionGC = gc_b, 
            percentile = percentiles, 
            transProb = keyObject$Transition_Probability, 
            chroms = NULL
        )
    })
    names(hmm) <- sampleIds
    hmm
}    
hmm <- reactive({
    req( isActiveVisibleTab(options) )
    batchEffects <- batchEffects()
    req(batchEffects)
    startSpinner(session, paste(module, 'getting HMM for each sample'))
    batchEffects$keyObject$Transition_Probability <- as.double( settings$Model_Parameters()$Transition_Probability$value ) # nolint
    hmm <- dataCache$get(CONSTANTS$compositeModelsType, keyObject = batchEffects$keyObject,
                         from = "ram", permanent = TRUE,
                         create = "asNeeded", createFn = createHMMs,
                         gcBias = gcBias(), batchEffects = batchEffects)
    stopSpinner(session)
    hmm    
})    

#----------------------------------------------------------------------
# create the key for cached plots
#----------------------------------------------------------------------
plotCacheKey <- reactive({
    list(
        batchEffectsKey = batchEffects()$key,
        settings = settings$Plot_Settings(),
        samples = sourceInfo()$key
    )
})

#----------------------------------------------------------------------
# plot results by genomic bin
#----------------------------------------------------------------------
getCNLines <- function(...){
    maxCN <- settings$Plot_Settings()$Maximum_Copy_Number$value
    color <- rep('grey', maxCN + 1)
    color[gcBiasSettings()$Fit_Parameters$Default_Modal_CN$value + 1] <- 'black'
    structure(0:maxCN, color = color)   
}
getChromLines <- function(d){
    L <- d$x[ copyNumber()[, isLastBin] ]
    structure(c(0, L) + 0.5, color = rep('grey', length(L) + 1))    
}
beCorrectedDepthPlotData <- reactive({
    copyNumber <- copyNumber()
    req(copyNumber)
    sampleIds <- clusteredSamples()
    d <- lapply(sampleIds, function(sampleId){
        data.table(
            x = copyNumber$binN,    
            y = copyNumber[[sampleId]],
            sampleId = sampleId
        )            
    })
    names(d) <- getSampleNames(sampleIds = sampleIds)
    d
})
hmmModelData <- reactive({
    hmm <- hmm()
    req(hmm)
    hmm <- hmm$value
    sampleIds <- clusteredSamples()
    d <- lapply(sampleIds, function(sampleId){
        cn_hmm <- hmm[[sampleId]]$cn
        nSegments <- length(cn_hmm$lengths)
        ends <- cumsum(cn_hmm$lengths) 
        starts <- c(1, ends + 1)[1:nSegments]        
        x <- sort(c(starts, ends))
        y <- as.vector(sapply(cn_hmm$values, rep, 2))
        data.table(
            x = x,
            y = y
        )         
    })
    names(d) <- getSampleNames(sampleIds = sampleIds)
    d
})
output$beCorrectedDepthPlotUI <- renderUI({
    interactiveScatterplotUI(ns('beCorrectedDepthPlot'),
                             height = paste0(length(goodSourceSamples()) *
                                           settings$Plot_Settings()$Track_Height_Pixels$value,
                                           "px"))
})
beCorrectedDepthPlot <- interactiveScatterplotServer(
    'beCorrectedDepthPlot',
    beCorrectedDepthPlotData,
    accelerate = TRUE,
    xtitle = "Genomic Bin",
    xrange = range_pad,
    xzeroline = FALSE,
    ytitle = "Copy Number",    
    yrange = reactive({ c(-0.25, settings$Plot_Settings()$Maximum_Copy_Number$value + 0.25) }),
    yzeroline = FALSE,
    shareAxis = list(x = TRUE),
    hLines = getCNLines,
    vLines = getChromLines,
    grid = list(x = FALSE, y = FALSE),
    cacheReactive = plotCacheKey,
    overplot = hmmModelData,
    overplotMode = 'lines',
    overplotColor = 'black',  
    overplotLineWidth = 2
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input, # inputs to be bookmarked as declared by boomarkedInput() wrapper  
    settings = settings$all_ # step settings for children and bookmarking
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
