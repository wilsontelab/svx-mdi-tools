
#----------------------------------------------------------------------
# support functions for normalizeGC
#----------------------------------------------------------------------

# recall or create the mdiSummarizedExperiment describing the input data
gcBiasIsReady <- function(keyObject, sampleIds=NULL){
    key <- dataCache$getCacheKeys(CONSTANTS$mseAssembledType, keyObject = keyObject)$key
    gcBiasOutcomes <- getStepOutcomesByType('normalizeGC')()
    req(gcBiasOutcomes[[key]])
    gcBias <- gcBiasOutcomes[[key]]$value
    libraryQCFailed <- getStepOutcomesByType('libraryQC')()
    batchReject <- gcBias$batchReject_s
    allSamplesIds <- names(batchReject)
    failed <- sapply(allSamplesIds, function(sampleId) libraryQCFailed[[keyObject$sourceId]][sampleId] )
    rejected <- as.integer(batchReject) == CONSTANTS$rejected     
    fitted   <- sapply(gcBias$fit, is.list) 
    required <- is.null(sampleIds) | allSamplesIds %in% sampleIds
    ready    <- failed | rejected | fitted | !required    
    sum(!ready) == 0 
}
getGcBiasMse <- function(module, options, sourceId, 
                         createFn=NULL, create="asNeeded",
                         sampleIds=NULL){
    req( isActiveVisibleTab(options) ) # only begin once the step/tab is in view
    req(sourceId)
    gcBiasSettings <- getStepSettingsByType('normalizeGC')()
    gcBiasSettings$Fit_Parameters$Default_Modal_CN <- NULL
    keyObject <- list(sourceId = sourceId, # define settings to pass, only the ones upon which we depend  
                      settings = list(Bin_Filters    = gcBiasSettings$Bin_Filters,
                                      Fit_Parameters = gcBiasSettings$Fit_Parameters))      
    req( !is.null(createFn) || gcBiasIsReady(keyObject, sampleIds) ) # for later app steps, demand that gcBias step was completed # nolint
    startSpinner(session, paste(module, 'getGcBiasMse'))
    mse <- dataCache$get(CONSTANTS$mseAssembledType, keyObject = keyObject, 
                         createFn = createFn, create = create,
                         from = "ram", permanent = TRUE)
    stopSpinner(session)
    mse
}

# recall the user outcomes for a specific keyed MSE
getGcBiasOutcomes <- function(mse){
    req(mse())
    gcBiasOutcomes <- getStepOutcomesByType('normalizeGC')()
    gcBias <- gcBiasOutcomes[[mse()$key]]
    req(gcBias)
    gcBias$value   
}

