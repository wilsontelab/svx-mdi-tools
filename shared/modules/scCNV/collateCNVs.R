#=====================================================================================
# compare CNVs across all cells in a project
#-------------------------------------------------------------------------------------
getCellCnvs <- function(shapeKey, cell){
    if(cell$badCell || is.null(cell$segments[[shapeKey]])) return(NULL)
    sampleName <- manifest[Sample_ID == cell$cell_id, Sample_Name]
    x <- copy(windows[[cell$windowPower + 1]])
    cellPloidy <- sampleChromosomes[[sampleName]]$CN[x$chrom, expected]
    x[, ":="(
        windowPower = cell$windowPower,
        keep = cell$keep,
        cellCN = cell$windows[[shapeKey]]$sequential$HMM,        
        sampleMedian = sampleProfiles[[sampleName]][[cell$windowPower + 1]],
        cellPloidy = cellPloidy
    )] 
    do.call(rbind, lapply(c("sampleMedian", "cellPloidy"), function(type){
        x[, .(
            type = type,
            sampleName = sampleName,
            cell_id = cell$cell_id, 
            chrom = chrom[1],
            start = min(start),
            end = max(end),
            gc_fraction = mean(gc_fraction),
            cellCN = cellCN[1],
            referenceCN = if(type == "sampleMedian") sampleMedian[1] else cellPloidy[1],
            nWindows = .N,
            windowPower = windowPower[1],
            keep = keep[1]
        ), by = rleidv(x, cols = c("chrom", "cellCN", type))][
            cellCN != referenceCN
        ]
    }))
}
collateCNVs <- function(){ # creates the unshaped and shaped, i.e., cell-specific, fits
    workingStep <<- "collateCNVs"
    shapeKeys <- c("unshaped", "shaped", "batched")
    cnvs <- lapply(shapeKeys, function(shapeKey){
        do.call(rbind, mclapply(cells, function(cell) {
            getCellCnvs(shapeKey, cell)
        }, mc.cores = env$N_CPU))
    })
    names(cnvs) <- shapeKeys
    cnvs        
}
#=====================================================================================
