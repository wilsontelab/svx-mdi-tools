#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# CONSTANTS
CONSTANTS$binSize <- 20000
CONSTANTS$nSdHalfCn <- 2
CONSTANTS$minMappability <- 0.9
CONSTANTS$minWindowPower <- 0
CONSTANTS$maxWindowPower <- 7
CONSTANTS$windowPowers <- CONSTANTS$minWindowPower:CONSTANTS$maxWindowPower
CONSTANTS$windowSizes <- 2 ** CONSTANTS$windowPowers
CONSTANTS$maxWindowSize <- 2 ** CONSTANTS$maxWindowPower

# cache for holding post-processed, but unchanging, cell-level data
windowsCache <- new_dataCache('scCNV-windows')
cellCache    <- new_dataCache('scCNV-cells')



# sample cache objects
sampleCache <- list(
    common = list( # rowRanges and colData, loaded once per sample
        cache = new_dataCache('scCNV-common'),
        claims = list()
    ), 
    working = list( # chrom-level expanded data, loaded per sample-chrom
        cache = new_dataCache('scCNV-working'),
        claims = list()
    )
)
getSampleCache <- function(caller, cacheType, sourceId, chrom = NULL){
    req(sourceId)
    dprint("getSampleCache")    
    key <- if(is.null(chrom)) sourceId else paste(sourceId, chrom, sep = "-")
    x <- sampleCache[[cacheType]]
    priorKeys <- x$cache$cacheKeys()
    d <- x$cache$get(
        cacheType, 
        key = key, 
        createFn = if(is.null(chrom)) loadSampleCommon else loadSampleWorking,
        permanent = FALSE,
        from = "ram",
        sourceId = sourceId,
        chrom = chrom
    )
    sampleCache[[cacheType]]$claims[[caller]] <<- d$cacheKey
    activeKeys <- unique(unlist(sampleCache[[cacheType]]$claims))
    staleKeys <- priorKeys[!(priorKeys %in% activeKeys)]
    for(key in staleKeys) x$cache$clear(key, purgeFiles = FALSE)
    d$value
}

# appServer function called after all modules are instantiated
appServer <- function(){


}
