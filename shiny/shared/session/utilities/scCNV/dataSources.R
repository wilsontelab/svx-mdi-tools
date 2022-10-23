projectNameReactive <- function(sourceId){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        gsub('.scCNV.extract', '', getSourceFilePackageName(sourceId))
    })    
} 
extractDataReactive <- function(sourceId){
    reactive({
        sourceId <- sourceId()
        req(sourceId)
        dataFilePath <- getSourceFilePath(sourceId, "extractFile")
        startSpinner(session, message = "loading sample data")
        x <- readRDS(dataFilePath)
        setkey(x$colData, "cell_id")  
        x$qcPlotsDir <- expandSourceFilePath(sourceId, "qc_plots")
        if(!dir.exists(x$qcPlotsDir)) untar(
            getSourceFilePath(sourceId, "plotsArchive"), 
            exdir = x$qcPlotsDir
        )
        stopSpinner(session)
        x
    }) 
}

# # expandWindows <- function(cell_id, data, key){ # go from window values to bin values, repeating window value over all bins
# #     cell <- cells[[cell_id]] # NOT data!
# #     w <- cell$window_size
# #     rollingRanges <- rollingRanges[[paste("w", w, sep = "_")]]
# #     chroms <- rollingRanges[reference_window == TRUE, chrom]
# #     unlist(rollingRanges[,
# #         sapply(data[[cell_id]][[key]][chroms == chrom[1]], function(x) rep(x, w))[1:.N],
# #         by = chrom
# #     ][, 2])
# # }
# # cn_w <- as.data.table(
# #     mclapply(constants$good_cell_ids, expandWindows, cells, "cn", mc.cores = env$N_CPU) 
# # )
