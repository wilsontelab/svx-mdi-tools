#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# objects instantiated here are available to all appStep modules in a session
svPoreCache <- new_dataCache('svPoreCache')
getJunctionClusters <- function(sourceId = NULL, item = NULL, coord = NULL, rangeType = NULL){
    if(is.null(sourceId)) {
        uploadName <- appStepNamesByType$upload
        samples <- as.data.table(app[[uploadName]]$outcomes$samples())        
        sourceId <- samples[Sample_ID == item$Sample_ID & Project == item$Project, Source_ID]
    }
    req(sourceId)
    jc <- svPoreCache$get('junctionClusters', key = sourceId, permanent = FALSE, from = "disk", create = "always", createFn = function(...) {
        startSpinner(session, message = "loading junction clusters")
        jc <- readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
        cm <- readRDS(getSourceFilePath(sourceId, "chromosomesFile")) 
        jc[, ":="(
            cChrom1 = unlist(cm$revChromIndex[cChromIndex1]),
            cChrom2 = unlist(cm$revChromIndex[cChromIndex2]),
            refPosCenter = pmin(cRefPos1, cRefPos2) + abs(cRefPos2 - cRefPos1) / 2,
            size = abs(cRefPos2 - cRefPos1)
        )]
        stopSpinner(session)
        jc
    })$value
    if(!is.null(coord) && coord$chromosome != "all"){
        jc <- jc[
            cChrom1 == coord$chromosome & 
            cChrom2 == coord$chromosome
        ]
    }
    if(!is.null(coord) && !is.null(rangeType)){
        jc <- switch(
            rangeType,
            center = jc[between(refPosCenter, coord$start, coord$end)],
            jc
        )
    }
    jc
}

# appServer function called after all modules are instantiated
appServer <- function(){


}
