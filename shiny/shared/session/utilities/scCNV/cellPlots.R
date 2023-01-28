#----------------------------------------------------------------------
# parse metadata on how to plot a cell given it's pipeline options and cell outcomes
#----------------------------------------------------------------------
shapeModels <- c( # in priority order
    "Batched",
    "Shaped",
    "Unshaped"
)
getShapeModel <- function(settings, cell){
    model <- settings$get("Page_Options", "Shape_Model")
    while(!(tolower(model) %in% names(cell$windows))) 
        model <- shapeModels[which(shapeModels == model) + 1]
    list(
        model = model,
        key = tolower(model)
    )
}
getReplicationModel <- function(settings, cell, getReplicating){
    model <- settings$get("Page_Options", "Replication_Model")
    isReplicating <- if(is.null(getReplicating)) cell$cellIsReplicating 
                     else getReplicating(cell$badCell, cell$cellIsReplicating, cell$cell_id)
    forceSequential <- isReplicating && model == "Sequential"    
    list(
        model = model,
        key = if(!isReplicating || forceSequential) "sequential" else "composite",
        isReplicating = isReplicating,
        forceSequential = forceSequential,
        isSequential = model == "Sequential"
    )
}

#----------------------------------------------------------------------
# make a composite plot of a single cell
#----------------------------------------------------------------------
pointOpacity <- function(x){
    min(1, max(0.15, -1/6000 * length(x) + 1))
}
getWindowColors <- function(x, pointOpacity = NULL){ # by CN (regardless of ploidy)
    if(is.null(pointOpacity)) pointOpacity <- pointOpacity(x)
    defaultPointColor <- rgb(0, 0, 0, pointOpacity)
    windowColors <- c(     
        defaultPointColor,                 # 0 = black/grey (absence of color/copies...)
        rgb(0.1, 0.1,  0.8, pointOpacity), # 1 = blue ("cool" colors are losses)
        rgb(0.1, 0.8,  0.1, pointOpacity), # 2 = subdued green ("good/go" for typical CN neutral)
        rgb(1,   0,    0,   pointOpacity), # 3 = red ("hot" colors are gains)
        rgb(1,   0.65, 0,   pointOpacity), # 4 = orange
        defaultPointColor                  # 5 = back to black/grey to make it obvious
    )
    col <- windowColors[pmin(x, 5) + 1]
    col[is.na(col)] <- defaultPointColor
    col
}
plotCellByWindow <- function(y, ylab, h, v, col = NULL, yaxt = NULL, ymax = NULL, cnvs = list()){
    if(is.null(col)) col <- rgb(0, 0, 0, pointOpacity(y))
    x <- 1:length(y)
    if(is.null(ymax)) ymax <- max(h)
    plot(NA, NA, bty = "n", xaxt = "n", yaxt = yaxt, xaxs = "i", , yaxs = "i",
         xlab = NULL, ylab = ylab, xlim = range(x), ylim = c(0, ymax))
    for(cnv in cnvs){
        rect(cnv$startI, 0, cnv$endI, ymax, col = getWindowColors(cnv$CN, 0.15), border = NA)
    }
    abline(h = h, v = v, col = "grey")
    points(x, y, pch = 19, cex = 0.4, col = col)
}
plotCellQC <- function(project, cell, shapeModel, replicationModel, settings){
    w <- project$windows[[cell$windowPower + 1]]
    v <- c(0, w[, max(i), by = "chrom"][[2]]) + 0.5
    if(cell$badCell){
        cw <- cell$windows$unshaped # bad cells only have a windows entry for windowPower 7, which is unshaped

        # plot NR_wms vs. gc_w
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        N <- length(cw$NR_wms)
        I <- sample.int(N, min(1e4, N))
        plot(w$gc_fraction[I], cw$NR_wms[I], 
            xlab = "Fraction GC", xlim = c(0.3, 0.6), 
            ylab = "# of Reads", #ylim = c(0, quantile(RPA, 0.975, na.rm = TRUE) * 1.5),
            pch = 19, cex = 0.4)

        # plot NR_wm vs. window index, i.e., pre-normalization input
        par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
        plotCellByWindow(cw$NR_wms, "# Reads", h = max(cw$NR_wms, na.rm = TRUE), v = v)
        
    } else {
        cw_pre <- cell$windows$unshaped
        cw_post <- cell$windows[[shapeModel$key]]
        cww <- cw_post[[replicationModel$key]]
        cww$NA_ <- cww$HMM + cww$NAR

        minWindows <- settings$get("Page_Options", "Minimum_Windows_Per_CNV", 20)
        cnvsRelativeTo <- settings$get("Page_Options", "CNVs_Relative_To", "Sample Median")
        cnvs <- project$cnvs[[shapeModel$key]][cell_id == cell$cell_id & nWindows >= minWindows]
        cnvs <- lapply(seq_len(nrow(cnvs)), function(j){
            x <- cnvs[j]
            w_i <- w[chrom == x$chrom & start ==x$start, i]
            dprint(paste(x$cellCN , cell$ploidy))
            list(
                startI = w_i,
                endI = w_i + x$nWindows - 1,
                CN = switch(
                    cnvsRelativeTo,
                    "Sample Median" = if(x$cellCN > x$sampleCN)  3 else 1,
                    "Cell Ploidy"   = if(x$cellCN > cell$ploidy) 3 else 1
                )
            )
        })

        # plot selected model's NR_wms vs. gc_w (color as final replication call)
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        N <- length(cw_post$NR_wms)
        N_d <- min(1e4, N)
        I <- sample.int(N, N_d)
        RPA <- cw_post$NR_wms[I] / (if(replicationModel$isSequential) cww$HMM[I] else cww$NA_[I])
        colNA <- if(replicationModel$isReplicating) round(cww$NAR[I] * cell$ploidy / cww$HMM[I], 0) else rep(0, N_d)
        plot(w$gc_fraction[I], RPA, 
            xlab = "Fraction GC", xlim = c(0.3, 0.6), 
            ylab = "Reads Per Allele", ylim = c(0, quantile(RPA[RPA < Inf], 0.975, na.rm = TRUE) * 1.5),
            pch = 19, cex = 0.4, 
            col = getWindowColors(colNA + 1)               
        )
        gc_fit <- if(replicationModel$isSequential) cww$gc_fit else cell$replicationModel[[shapeModel$key]]$gc_fit 
        with(gc_fit, { for(percentile in c(0.025, 0.975)) lines(
            gcFractions, 
            qnbinom(percentile, size = theta, mu = mu), 
            lty = 3, lwd = 1.5, col = "grey30"
        ) })
        
        # plot NR_wm vs. window index, i.e., pre-normalization input (color as final replication call)
        par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
        maxPlotNa <- max(cww$NA_, na.rm = TRUE) + 2
        colNA <- if(replicationModel$isReplicating) round(cww$NAR * cell$ploidy / cww$HMM, 0) else rep(0, length(cww$HMM))
        plotCellByWindow(cw_pre$NR_wms, "# Reads", 
                         h = cw_pre$RPA * 0:round(maxPlotNa, 0), v = v, 
                         ymax = max(cw_pre$RPA * maxPlotNa, quantile(cw_pre$NR_wms, 0.95, na.rm = TRUE)),
                         col = getWindowColors(colNA + 1), cnvs = cnvs)

        # plot CN vs. window index, i.e., post-normalization (color as final CN call)
        par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
        maxPlotNa <-  max(cww$HMM, 3, na.rm = TRUE) + 1
        sampleName <- project$manifest[Sample_ID == cell$cell_id, Sample_Name]
        sampleHMM <- project$sampleProfiles[[sampleName]][[cell$windowPower + 1]]
        plotCellByWindow(cww$CN, "CN", h = 0:maxPlotNa, v = v, 
                         col = getWindowColors(cww$HMM - sampleHMM + cell$ploidy), 
                         yaxt = "n", cnvs = cnvs)
        axis(2, at = 0:maxPlotNa, labels = 0:maxPlotNa) 

        # add HMM traces to CN plot
        lines(1:length(sampleHMM), sampleHMM, col = "black")
        segments <- cell$segments[[shapeModel$key]][[replicationModel$key]]$CN
        maskedHMM <- unlist(sapply(1:nrow(segments), function(j){
            x <- segments[j]
            length <- x[, chromEndI - chromStartI + 1]
            c_i <- x$chromStartI:x$chromEndI
            w_i <- w[chrom == x$chrom, i]
            cell   <- cww$HMM[w_i][c_i]
            sample <- sampleHMM[w_i][c_i]

            tooSmall <- length < minWindows # TODO:1) cnv areas 2) save outcomes,etc. 

            ifelse(cell == sample, cell, if(tooSmall) rep(as.integer(NA), length) else cell)
        }))
        lines(1:length(maskedHMM), maskedHMM, col = "brown")
    }
}
getCellCompositePlot <- function(project, cell, settings, getReplicating){
    name <- if(cell$badCell) "bad" else {
        shapeModel <- getShapeModel(settings, cell)
        replicationModel <- getReplicationModel(settings, cell, getReplicating)
        paste(shapeModel$key, replicationModel$key, replicationModel$isReplicating, sep = "_")
    }
    fileName <- paste(cell$cell_id, "qc", name, "png", sep = ".")
    pngFile <- file.path(project$qcPlotsDir, fileName)

    #####################
    # unlink(file.path(project$qcPlotsDir, "*.png"))
    force <- TRUE

    if(force || !file.exists(pngFile)) {
        png(
            filename = pngFile,
            width  = 1.5 * 6, 
            height = 1.35, 
            units = "in", 
            pointsize = 7,
            bg = "white",
            res = 96, # i.e., optimized for on-screen display
            type = "cairo"
        )
        layout(matrix(c(c(1,1), rep(c(2,3), 5)), nrow = 2, ncol = 6))
        plotCellQC(project, cell, shapeModel, replicationModel, settings)
        dev.off()
    }
    tags$img(src = pngFileToBase64(pngFile), style = "vertical-align: top;")
}

#----------------------------------------------------------------------
# create a div of metadata about a cell and buttons to manipulate it
#----------------------------------------------------------------------
getCellSummary <- function(project, cell, buttons){
    colData <- project$colData[cell$cell_id]
    desc <- project$manifest[as.character(Sample_ID) == cell$cell_id, Description]
    tags$div(
        style = "display: inline-block; padding: 5px; position: relative; z-index: 999; height: 130px; background: white;",
        tags$div(tags$strong( paste("cell:", colData$cell_id, '(', desc, ')') )), 
        tags$div(paste("windowPower: ", colData$windowPower, '(', commify(2 ** colData$windowPower * 20), 'kb )')),
        if(cell$badCell) tags$div("bad cell, not processed") else tagList(
            tags$div(paste("repGcRatio:", round(colData$repGcRatio, 2), '(', colData$modelType, ',', round(colData$fractionS, 2), ')')), 
            tags$div(paste("cnsd:", round(colData$cnsd, 2), '(', if(colData$keep) "" else "not", "passed", ')')),
            if(!is.null(buttons)) buttons(colData, cell) else "" 
        ) 
    )
}

#----------------------------------------------------------------------
# master function for plotting a single cell
#----------------------------------------------------------------------
plotOneCellUI <- function(project, cell, settings, buttons = NULL, getKeep = NULL, getReplicating = NULL){
    req(project)
    col <- if(cell$badCell || is.null(getKeep)) "grey" 
           else if(getKeep(cell$badCell, cell$keep, cell$cell_id)) "rgb(0,0,175)" 
           else "rgb(150,0,0)"
    tags$div(
        style = paste("border: 1px solid", col, "; white-space: nowrap;"),
        tagList(
            getCellCompositePlot(project, cell, settings, getReplicating),
            getCellSummary(project, cell, buttons)                
        )
    )
}
