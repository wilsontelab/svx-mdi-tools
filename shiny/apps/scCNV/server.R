#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# help manage the server (not session) scCNV data cache
getScCnvProjectData <- function(sourceId){
    startSpinner(session, message = "loading project data")
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "normalizeFile", 
        ttl = 24 * 60 * 60,
        postProcess = function(x){
            setkey(x$colData, "cell_id")  
            x$qcPlotsDir <- expandSourceFilePath(sourceId, "qc_plots")
            if(!dir.exists(x$qcPlotsDir)) dir.create(x$qcPlotsDir)
            x$rowRanges <- NULL # remove the objects that inflate the size of the full project file for faster caching
            x$raw_counts <- NULL
            x
        }
    )
    stopSpinner(session)
    persistentCache[[filePath]]$data
}
scCnvProjectReactive <- function(sourceId, ...){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        getScCnvProjectData(sourceId)
    }) 
}

# # sample cache objects
# sampleCache <- list(
#     common = list( # rowRanges and colData, loaded once per sample
#         cache = new_dataCache('scCNV-common'),
#         claims = list()
#     ), 
#     working = list( # chrom-level expanded data, loaded per sample-chrom
#         cache = new_dataCache('scCNV-working'),
#         claims = list()
#     )
# )

# # cache for holding post-processed, but unchanging, cell-level data
# windowsCache <- new_dataCache('scCNV-windows')
# cellCache    <- new_dataCache('scCNV-cells')

# # sample cache objects
# sampleCache <- list(
#     common = list( # rowRanges and colData, loaded once per sample
#         cache = new_dataCache('scCNV-common'),
#         claims = list()
#     ), 
#     working = list( # chrom-level expanded data, loaded per sample-chrom
#         cache = new_dataCache('scCNV-working'),
#         claims = list()
#     )
# )
# getSampleCache <- function(caller, cacheType, sourceId, chrom = NULL){
#     req(sourceId)
#     dprint("getSampleCache")    
#     key <- if(is.null(chrom)) sourceId else paste(sourceId, chrom, sep = "-")
#     x <- sampleCache[[cacheType]]
#     priorKeys <- x$cache$cacheKeys()
#     d <- x$cache$get(
#         cacheType, 
#         key = key, 
#         createFn = if(is.null(chrom)) loadSampleCommon else loadSampleWorking,
#         permanent = FALSE,
#         from = "ram",
#         sourceId = sourceId,
#         chrom = chrom
#     )
#     sampleCache[[cacheType]]$claims[[caller]] <<- d$cacheKey
#     activeKeys <- unique(unlist(sampleCache[[cacheType]]$claims))
#     staleKeys <- priorKeys[!(priorKeys %in% activeKeys)]
#     for(key in staleKeys) x$cache$clear(key, purgeFiles = FALSE)
#     d$value
# }

# appServer function called after all modules are instantiated
appServer <- function(){


}
