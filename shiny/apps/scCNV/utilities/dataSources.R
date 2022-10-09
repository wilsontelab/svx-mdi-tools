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
        startSpinner(session, message = "loading sample data")
        x <- readRDS(dataFilePath)
        setkey(x$colData, "cell_id")
        x$chromEnds <- x$rowRanges[, max(bin_n, na.rm = TRUE), by = chrom][[2]]
        stopSpinner(session)
        x
    }) 
}
