#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# create constants that define use actions on cells
CONSTANTS$batched   <- 0 # cell/sample used to determine batch effects
CONSTANTS$unbatched <- 1 # cell/sample omitted from batch effect calculations but shown in composites
CONSTANTS$rejected  <- 2 # cell/sample never used after gcBias correction
CONSTANTS$noCnCall  <- -99 # a flag that a bin or run has not yet been assigned a CN state
getSampleFlagStatus <- Vectorize(function(flag){
    flag <- as.integer(flag)
    if(flag == CONSTANTS$noCnCall) '_pending_'
    else switch(flag + 1, 'Batch', 'No Batch', 'Rejected')
})

# initialize data cache, available to all modules in this session
dataCache <- new_dataCache('SampleSets')
CONSTANTS$mseInitialType     <- 'mdiSummarizedExperiment_initial'
CONSTANTS$mseFilteredType    <- 'mdiSummarizedExperiment_filtered'
CONSTANTS$mseAssembledType   <- 'mdiSummarizedExperiment_asembled'
CONSTANTS$gcBiasOutcomesType <- 'gcBiasOutcomes'
CONSTANTS$batchEffectOutcomesType   <- 'batchEffectOutcomes'
CONSTANTS$compositeModelsType       <- 'compositeModels'
CONSTANTS$compositeCopyNumberType   <- 'compositeCopyNumber'

# appServer function called after all modules are instantiated
appServer <- function(){

}
