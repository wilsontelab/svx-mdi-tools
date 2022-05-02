#----------------------------------------------------------------------
# support functions for normalizeBatch
#----------------------------------------------------------------------
getBatchEffects <- function(module, options,
                            mse, gcBias, sampleIds, 
                            createFn=NULL, create="asNeeded"){
    req( isActiveVisibleTab(options) )
    req(mse$value)
    req(gcBias)
    startSpinner(session, paste(module, 'getting batchEffects'))    
    keyObject <- mse$keyObject # static and gcBias outcomes objects are keyed the same
    keyObject$gcBias <- gcBias[c('batchReject_s', 'MCN_s', 'hmm')] # as well as the specific outcomes of the gcBias step
    keyObject$sampleIds <- sort(unique(sampleIds)) # and the subset of samples in this batch
    be <- dataCache$get(CONSTANTS$batchEffectOutcomesType, keyObject = keyObject,
                        from = "ram", permanent = TRUE,
                        create = create, createFn = createFn,
                        mse = mse, gcBias = gcBias, sampleIds = sampleIds)
    stopSpinner(session)
    req(be)
    be 
}
