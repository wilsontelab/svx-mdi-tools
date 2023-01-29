#----------------------------------------------------------------------
# parse settings on how to plot a cell given it's pipeline options and cell outcomes
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
getCnvsType <- function(settings){
    cnvsType <- settings$get("Page_Options", "CNVs_Relative_To", "Sample Median")
    if(cnvsType == "Sample Median") "sampleMedian" else "cellPloidy" 
}

#----------------------------------------------------------------------
# window point colors
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
adjustColorAlpha <- function(col, minAlpha) {
    alpha <- max(minAlpha, pointOpacity(col))
    apply(sapply(col, col2rgb) / 255, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}
#----------------------------------------------------------------------
# make a composite plot of a single cell (on a single chromosome)
#----------------------------------------------------------------------
plotCellByWindow <- function(y, ylab, h, v, col = NULL, yaxt = NULL, ymax = NULL, 
                             cnvs = list(), sampleHMM = NULL, maskedHMM = NULL, chromI = NULL){
    if(is.null(ymax)) ymax <- max(h)
    yInset <- 0 # ymax * 0.01 
    isZoomed <- !is.null(chromI)
    I <- if(isZoomed) chromI else TRUE
    x <- 1:length(y) 
    if(isZoomed){
        x <- x[I]
        y <- y[I]
    }
    if(is.null(col)) col <- rgb(0, 0, 0, pointOpacity(y))
    else if(isZoomed) col <- adjustColorAlpha(col[I], 0.5)
    par(mar = c(0.1, if(isZoomed) 0 else 4.1, 0.1, if(isZoomed) 0 else 0.1), cex = 1)  
    plot(NA, NA, bty = "n", xaxt = "n", yaxt = if(isZoomed) "n" else yaxt, xaxs = "i", , yaxs = "i",
        xlab = NULL, ylab = ylab, xlim = range(x), ylim = c(0, ymax))
    for(cnv in cnvs){
        rect(
            cnv$startI - 0.5, yInset, cnv$endI + 0.5, ymax - yInset, 
            col = getWindowColors(cnv$CN, 0.15), 
            border = NA
        )
    }
    abline(h = h, v = v, col = "grey")
    points(x, y, pch = 19, cex = 0.4, col = col)
    if(!is.null(sampleHMM)) lines(x, sampleHMM[I], col = "black")        
    if(!is.null(maskedHMM)) lines(x, maskedHMM[I], col = "brown")
    if(!isZoomed && !is.null(yaxt)) axis(2, at = 0:ymax, labels = 0:ymax) 
}
plotCellQC <- function(cell, d, zoomChrom){
    isZoomed <- !is.null(zoomChrom)
    zc_i <- if(isZoomed) d$w[chrom == zoomChrom, i] else NULL
    if(cell$badCell){
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        if(!isZoomed) plot(d$w$gc_fraction[d$d$I_d], d$d$cw$NR_wms[d$d$I_d], 
            xlab = "Fraction GC", xlim = c(0.3, 0.6), 
            ylab = "# of Reads", 
            pch = 19, cex = 0.4)
        plotCellByWindow(d$d$cw$NR_wms, "# Reads", h = max(d$d$cw$NR_wms, na.rm = TRUE), v = d$v, chromI = zc_i)
    } else {

        # plot selected model's NR_wms vs. gc_w (color as final replication call)
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        if(!isZoomed) {
            plot(d$w$gc_fraction[d$d$I_d], d$d$RPA, 
                xlab = "Fraction GC", xlim = c(0.3, 0.6), 
                ylab = "Reads Per Allele", ylim = c(0, quantile(d$d$RPA[d$d$RPA < Inf], 0.975, na.rm = TRUE) * 1.5),
                pch = 19, cex = 0.4, 
                col = d$d$colNA_gc          
            )
            gc_fit <- if(d$replicationModel$isSequential) d$d$cww$gc_fit else cell$replicationModel[[d$shapeModel$key]]$gc_fit 
            with(gc_fit, { for(percentile in c(0.025, 0.975)) lines(
                gcFractions, 
                qnbinom(percentile, size = theta, mu = mu), 
                lty = 3, lwd = 1.5, col = "grey30"
            ) })
        }        

        # plot NR_wm vs. window index, i.e., pre-normalization input (color as final replication call) 
        plotCellByWindow(d$d$cw_pre$NR_wms, "# Reads", h = d$d$h_NR, v = d$v, ymax = d$d$ymax_NR,
                         col = d$d$colNA_NR, cnvs = d$d$cnvs, chromI = zc_i)

        # plot CN vs. window index, i.e., post-normalization (color as final CN call)
        plotCellByWindow(d$d$cww$CN, "CN", h = d$d$h_CN, v = d$v, yaxt = "n", 
                         col = d$d$colNA_CN, cnvs = d$d$cnvs, chromI = zc_i, 
                         sampleHMM = d$d$sampleHMM, maskedHMM = d$d$maskedHMM)
    }
}
getCellCompositeKey <- function(shapeModel, replicationModel, settings){
    paste(
        shapeModel$key, 
        replicationModel$key, 
        replicationModel$isReplicating,
        settings$get("Page_Options", "Minimum_Windows_Per_CNV", 20),
        getCnvsType(settings),
        sep = "_"
    )
}
getCellPngFile <- function(project, cell, level, d){
    fileName <- paste(cell$cell_id, level, d$key, "png", sep = ".")
    file.path(project$qcPlotsDir, fileName)
}
getCellPlotPanel <- function(project, cell, settings, level, 
                             widthChunks, layoutMatrix, zoomChrom = NULL, force = TRUE){
    d <- cellCompositeData[[cell$cell_id]]
    pngFile <- getCellPngFile(project, cell, level, d)
    mustCreate <- force || !file.exists(pngFile)
    if(mustCreate) {
        png(
            filename = pngFile,
            width  = 1.5 * widthChunks, 
            height = 1.35, 
            units = "in", 
            pointsize = 7,
            bg = "white",
            res = 96, # i.e., optimized for on-screen display
            type = "cairo"
        )
        layout(layoutMatrix)
        plotCellQC(cell, d, zoomChrom)
        dev.off()
    }
    tags$img(
        src = pngFileToBase64(pngFile), 
        style = "vertical-align: top; cursor: pointer;", 
        class = "cellWindowsPlot",
        "data-cell_id" = cell$cell_id,
        "data-level" = if(level == "genome") "genome" else "chrom"
    )
}
getCellCompositePlot <- function(project, cell, settings){

    #####################
    force <- FALSE    

    getCellPlotPanel(
        project, cell, settings, "genome", 
        widthChunks = 6,
        layoutMatrix = matrix(c(c(1,1), rep(c(2,3), 5)), nrow = 2, ncol = 6),
        force = force
    )  
}
getZoomChromPlot <- function(project, cell, settings, zoomChrom){
    if(is.null(zoomChrom)) return("")
    getCellPlotPanel(
        project, cell, settings, zoomChrom, 
        widthChunks = 3,
        layoutMatrix = matrix(rep(c(1,2), 3), nrow = 2, ncol = 3),
        zoomChrom = zoomChrom,
        force = FALSE
    )  
}
setCellCompositeData <- function(project, cell, settings, level,
                                 getReplicating = NULL, getKeep = NULL){
    key <- if(cell$badCell){
        shapeModel <- NULL
        replicationModel <- NULL
        "badCell"   
    } else {
        shapeModel <- getShapeModel(settings, cell)
        replicationModel <- getReplicationModel(settings, cell, getReplicating)
        getCellCompositeKey(shapeModel, replicationModel, settings)
    }
    d <- list(key = key)

    divCol <- if(cell$badCell || is.null(getKeep)) "grey" 
         else if(getKeep(cell$badCell, cell$keep, cell$cell_id)) "rgb(0,0,175)" 
         else "rgb(150,0,0)"

    pngFile <- getCellPngFile(project, cell, level, d)
    imageFileExists <- file.exists(pngFile)
    keyIsLoaded <- !is.null(cellCompositeData[[cell$cell_id]]) &&
                   cellCompositeData[[cell$cell_id]]$key == key
    dataIsLoaded <- keyIsLoaded && !is.null(cellCompositeData[[cell$cell_id]]$d)
    if(imageFileExists){
        if(!keyIsLoaded) cellCompositeData[[cell$cell_id]] <<- list(
            key = key,
            shapeModel = shapeModel,
            replicationModel = replicationModel,
            col = divCol
        )
        return(NULL)
    } else if(dataIsLoaded){
        return(NULL)
    }

    w <- project$windows[[cell$windowPower + 1]]
    v <- c(0, w[, max(i), by = "chrom"][[2]]) + 0.5
    d <- if(cell$badCell){
        cw <- cell$windows$unshaped # bad cells only have a windows entry for windowPower 7, which is unshaped
        N <- length(cw$NR_wms)
        list(
            cw = cw,
            I_d = sample.int(N, min(1e4, N))
        )
    } else {
        cw_pre <- cell$windows$unshaped
        cw_post <- cell$windows[[shapeModel$key]]
        cww <- cw_post[[replicationModel$key]]
        cww$NA_ <- cww$HMM + cww$NAR

        N <- length(cw_post$NR_wms)
        N_d <- min(1e4, N)
        I_d <- sample.int(N, N_d)  

        maxPlotNa_NR <- max(cww$NA_, na.rm = TRUE) + 2
        maxPlotNa_CN <-  max(cww$HMM, 3, na.rm = TRUE) + 1

        RPA <- cw_post$NR_wms[I_d] / (if(replicationModel$isSequential) cww$HMM[I_d] else cww$NA_[I_d])
        gc_fit <- if(replicationModel$isSequential) cww$gc_fit else cell$replicationModel[[shapeModel$key]]$gc_fit

        minWindows <- settings$get("Page_Options", "Minimum_Windows_Per_CNV", 20)
        cnvsType <- getCnvsType(settings)
        cnvs <- project$cnvs[[shapeModel$key]][type == cnvsType & cell_id == cell$cell_id & nWindows >= minWindows]
        cnvs <- lapply(seq_len(nrow(cnvs)), function(j){
            x <- cnvs[j]
            w_i <- w[chrom == x$chrom & start ==x$start, i]
            list(
                startI = w_i,
                endI = w_i + x$nWindows - 1,
                CN = if(x$cellCN > x$referenceCN) 3 else 1
            )
        })

        segments <- cell$segments[[shapeModel$key]][[replicationModel$key]]$CN        
        sampleName <- project$manifest[Sample_ID == cell$cell_id, Sample_Name]
        sampleHMM <- project$sampleProfiles[[sampleName]][[cell$windowPower + 1]]
        referenceCN <- if(cnvsType == "sampleMedian") sampleHMM else cell$ploidy        
        maskedHMM <- unlist(sapply(1:nrow(segments), function(j){
            x <- segments[j]
            length <- x[, chromEndI - chromStartI + 1]
            c_i <- x$chromStartI:x$chromEndI
            w_i <- w[chrom == x$chrom, i]
            cellCN   <- cww$HMM[w_i][c_i]
            referenceCN <- if(cnvsType == "sampleMedian") sampleHMM[w_i][c_i] else cell$ploidy
            tooSmall <- length < minWindows
            ifelse(cellCN == referenceCN, cellCN, if(tooSmall) rep(as.integer(NA), length) else cellCN)
        }))
        list(
            cw_pre = cw_pre,
            cw_post = cw_post,
            cww = cww,
            N_d = N_d,
            I_d = I_d,
            colNA_gc = getWindowColors((if(replicationModel$isReplicating) round(cww$NAR[I_d] * cell$ploidy / cww$HMM[I_d], 0) 
                                       else rep(0, N_d)) + 1),
            colNA_NR = getWindowColors((if(replicationModel$isReplicating) round(cww$NAR * cell$ploidy / cww$HMM, 0) 
                                       else rep(0, length(cww$HMM))) + 1),
            colNA_CN = getWindowColors(cww$HMM - referenceCN + cell$ploidy),
            maxPlotNa_NR = maxPlotNa_NR,
            maxPlotNa_CN = maxPlotNa_CN,
            ymax_NR = max(cw_pre$RPA * maxPlotNa_NR, quantile(cw_pre$NR_wms, 0.95, na.rm = TRUE)),
            h_NR = cw_pre$RPA * 0:round(maxPlotNa_NR, 0),
            h_CN = 0:maxPlotNa_CN, 
            RPA = RPA,
            gc_fit = gc_fit,
            cnvs = cnvs,
            referenceCN = referenceCN,
            ploidy = cell$ploidy,
            sampleHMM = sampleHMM,
            maskedHMM = maskedHMM
        )
    }
    cellCompositeData[[cell$cell_id]] <<- list(
        key = key,
        shapeModel = shapeModel,
        replicationModel = replicationModel,
        col = divCol,
        w = w,
        v = v,
        d = d
    )
}

#----------------------------------------------------------------------
# create a div of metadata about a cell and buttons to manipulate it
#----------------------------------------------------------------------
getCellSummary <- function(project, cell, buttons){
    colData <- project$colData[cell$cell_id]
    desc <- project$manifest[as.character(Sample_ID) == cell$cell_id, Description]
    tags$div(
        class = "cellSummaryPanel",
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
cellCompositeData <- list()
plotOneCellUI_genome <- function(project, cell, settings, buttons = NULL, 
                                 getReplicating = NULL, getKeep = NULL){
    req(project)
    unlink(file.path(project$qcPlotsDir, "*.png"))
    setCellCompositeData(project, cell, settings, "genome", getReplicating, getKeep)
    tags$div(
        style = paste("border: 1px solid", cellCompositeData[[cell$cell_id]]$col, ";"),
        class = "cellPlotWrapper",
        tagList(
            getCellSummary(project, cell, buttons),          
            getCellCompositePlot(project, cell, settings)
        )
    )
}
plotOneCellUI_chrom <- function(project, cell, settings, zoomChrom,
                                 getReplicating = NULL, getKeep = NULL){
    req(project)
    setCellCompositeData(project, cell, settings, zoomChrom, getReplicating, getKeep)
    tags$div(
        style = paste("border: 1px solid", cellCompositeData[[cell$cell_id]]$col, ";"),
        class = "cellPlotWrapper",
        getZoomChromPlot(project, cell, settings, zoomChrom)
    )
}

#----------------------------------------------------------------------
# handle plot click
#----------------------------------------------------------------------
genomeXOffset <- (1.5 + 0.5) * 96
genomeXWidth <- (1.5 * 5 - 0.5) * 96
chromXOffset <- 0
chromXWidth <- (1.5 * 3 - 0) * 96
handleCellPlotClick <- function(project, click, zoomChrom, zoomWindowI){
    expandClick <- function(offest, width){
        x <- click$coord$x - offest
        cell_id <- as.character(click$data$cell_id) # gets sent to us as an integer by Shiny.setInputValue
        cell <- project$cells[[cell_id]]
        w <- project$windows[[cell$windowPower + 1]] 
        list(
            x = x,
            w = w
        )       
    }
    switch(
        click$data$level,
        genome = {
            click <- expandClick(genomeXOffset, genomeXWidth)
            chrom <- click$w[round(click$x / genomeXWidth * nrow(click$w), 0), chrom]
            zoomChrom(chrom)
        },
        chrom = {
            click <- expandClick(chromXOffset, chromXWidth)
            w_c <- click$w[chrom == zoomChrom()]
            n_w <- nrow(w_c)
            I <- round(click$x / chromXWidth * n_w, 0)
            I <- pmin(n_w, pmax(1, I))
            zoomWindowI(w_c[I, i])
        }
    )
}
