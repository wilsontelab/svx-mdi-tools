#=====================================================================================
# main target function that compares CNVs across all cells in a project
#-------------------------------------------------------------------------------------
collateCNVs <- function(cells){ # creates the unshaped and shaped, i.e., cell-specific, fits
    # plotNumber <<- 1 # so developer plots are ordered in lists as they were created
    workingStep <<- "collateCNVs"
    shapeKeys <- c("unshaped", "shaped", "batched")
    shapeKeys <- shapeKeys[shapeKeys %in% names(cells[[1]]$segments)]
    lapply(shapeKeys, function(shapeKey){
        do.call(rbind, mclapply(cells, function(cell) {
            if(cell$badCell) return(NULL)
            sampleName <- manifest[Sample_ID == cell$cell_id, Sample_Name]
            x <- copy(windows[[cell$windowPower + 1]])
            x[, ":="(
                windowPower = cell$windowPower,
                keep = cell$keep,
                sampleCN = sampleProfiles[[sampleName]][[cell$windowPower + 1]],
                cellCN = cell$windows[[shapeKey]]$sequential$HMM
            )]
            x[, .(
                sampleName = sampleName,
                cell_id = cell$cell_id, 
                chrom = chrom[1],
                start = min(start),
                end = max(end),
                gc_fraction = mean(gc_fraction),
                sampleCN = sampleCN[1],
                cellCN = cellCN[1],
                nWindows = .N,
                windowPower = windowPower[1],
                keep = keep[1]
            ), by = rleidv(x, cols = c("chrom", "sampleCN", "cellCN"))][sampleCN != cellCN]
        }, mc.cores = env$N_CPU))
    })
}
#=====================================================================================
