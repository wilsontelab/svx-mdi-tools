projectNameReactive <- function(sourceId){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        gsub('.scCNV.do', '', getSourceFilePackageName(sourceId))
    })    
}
sampleDataReactive <- function(sourceId){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        dataFilePath <- getSourceFilePath(sourceId, "normalizeFile")
        readRDS(dataFilePath)
    }) 
}
