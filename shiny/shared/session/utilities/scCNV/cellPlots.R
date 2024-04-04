#----------------------------------------------------------------------
# parse settings on how to plot a cell given it's pipeline options and cell outcomes
#----------------------------------------------------------------------
shapeModels <- c( # in priority order
    "Batched",
    "Shaped",
    "Unshaped"
)
getShapeModel <- function(settings, cell = NULL, cells = NULL){
    if(is.null(cell)){
        goodCellI <- 1
        while(cells[[goodCellI]]$badCell) goodCellI <- goodCellI + 1    
        cell <- cells[[goodCellI]]    
    }
    model <- settings$get("Page_Options", "Shape_Model")
    while(!(tolower(model) %in% names(cell$windows))) 
        model <- shapeModels[which(shapeModels == model) + 1]
    list(
        model = model,
        key = tolower(model)
    )
}
getReplicationModel <- function(settings, cell, getReplicating = NULL){
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
getCnvKey <- function(cnv, sourceId){ # a single row of a cnvs data.table
    cnv[, paste(type, sourceId, cell_id, chrom, start)]
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
plotCellByWindow <- function(y, ylab, h, v, isKeptCnv, isMatchedCnv, col = NULL, yaxt = NULL, ymax = NULL, 
                             cnvs = list(), sampleHMM = NULL, maskedHMM = NULL, chromI = NULL){
    if(is.null(ymax)) ymax <- max(h)
    yInset <- ymax * 0.025
    isZoomed <- !is.null(chromI)
    areMatching <- !is.null(isMatchedCnv)
    I <- if(isZoomed) chromI else TRUE
    x <- 1:length(y) 
    if(isZoomed){
        x <- x[I]
        y <- y[I]
    }
    if(is.null(col)) col <- rgb(0, 0, 0, pointOpacity(y))
    else if(isZoomed) col <- adjustColorAlpha(col[I], 0.5)
    par(mar = c(0.1, if(isZoomed) 0 else 4.1, 0.1, if(isZoomed) 0 else 0.1), cex = 1)  
    plot(NA, NA, bty = "n", xaxt = "n", yaxt = if(isZoomed) "n" else yaxt, xaxs = "i", yaxs = "i",
        xlab = NULL, ylab = ylab, xlim = range(x), ylim = c(0, ymax))

    addCnvRect <- function(cnv){
        # isColoredCnv <- !isZoomed || if(areMatching) isMatchedCnv(key = cnv$key) else isKeptCnv(key = cnv$key)
        cnvMatch <- if(areMatching) isMatchedCnv(cnvKey = cnv$key) else NULL
        isColoredCnv <- !isZoomed || if(areMatching) NA else isKeptCnv(key = cnv$key)
        cnvCol <- if(areMatching) getWindowColors(cnvMatch$CN_col, 0.15)
                  else if(isColoredCnv) getWindowColors(cnv$CN, 0.15)
                  else getWindowColors(0, 0.15)
        borderCol <- if(areMatching) { if(cnvMatch$isAssigned) getWindowColors(cnvMatch$CN_border, 0.8) else NA }
                     else if(isZoomed && isColoredCnv) getWindowColors(cnv$CN, 0.8)
                     else NA
        rect(
            cnv$startI, yInset, cnv$endI, ymax - yInset, 
            col = cnvCol, 
            border = borderCol,
            lwd = 2,
            lty = if(areMatching) cnvMatch$lty else 1
        )
    }
    if(is.data.table(cnvs)) for(i in seq_len(nrow(cnvs))) addCnvRect(cnvs[i, ])
    else for(cnv in cnvs) addCnvRect(cnv)
    abline(h = h, v = v, col = "grey")
    points(x, y, pch = 19, cex = 0.4, col = col)
    if(!is.null(sampleHMM)) lines(x, sampleHMM[I], col = "black")        
    if(!is.null(maskedHMM)) lines(x, maskedHMM[I], col = "brown")
    if(!isZoomed && !is.null(yaxt)) axis(2, at = 0:ymax, labels = 0:ymax) 
}
plotCellQC <- function(cell, d, zoomChrom, isKeptCnv, isMatchedCnv, cnvs, short){
    isZoomed <- !is.null(zoomChrom)
    zc_i <- if(isZoomed) d$w[chrom == zoomChrom, i] else NULL
    if(cell$badCell){
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        if(!isZoomed) plot(d$w$gc_fraction[d$d$I_d], d$d$cw$NR_wms[d$d$I_d], 
            xlab = "Fraction GC", xlim = c(0.3, 0.6), 
            ylab = "# of Reads", 
            pclickch = 19, cex = 0.4)
        plotCellByWindow(d$d$cw$NR_wms, "# Reads", h = max(d$d$cw$NR_wms, na.rm = TRUE), v = d$v, 
                         isKeptCnv = isKeptCnv, isMatchedCnv = isMatchedCnv, chromI = zc_i)
    } else {

        # plot selected model's NR_wms vs. gc_w (color as final replication call)
        par(mar = c(4.1, 4.1, 0.1, 1.1), cex = 1)
        if(!isZoomed && !short) {
            plot(d$w$gc_fraction[d$d$I_d], d$d$RPA_post_d, 
                xlab = "Fraction GC", xlim = c(0.3, 0.6), 
                ylab = "Reads Per Allele", ylim = c(0, quantile(d$d$RPA_post_d[d$d$RPA_post_d < Inf], 0.975, na.rm = TRUE) * 1.5),
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
        if(is.null(cnvs)) cnvs <- d$d$cnvs
        # if(!short) {
            sexChromModeHmm <- mode(d$d$maskedHMM[!d$w$autosome])
            expectedCn <- ifelse(d$w$autosome, cell$ploidy, sexChromModeHmm)
            nrCol <- if(cell$cellIsReplicating) d$d$colNA_NR else ifelse(
                is.na(d$d$cww$HMM) | d$d$cww$HMM == expectedCn, 
                rgb(0.2, 0.2, 0.2, pointOpacity(expectedCn)), 
                d$d$colNA_CN
            )
            # plotCellByWindow(d$d$cw_pre$NR_wms, "# Reads", h = d$d$h_NR, v = d$v, 
            #              isKeptCnv = isKeptCnv, isMatchedCnv = isMatchedCnv, ymax = d$d$ymax_NR,
            #              col = nrCol, cnvs = cnvs, chromI = zc_i)
            plotCellByWindow(d$d$cn_pre_w, "NR / RPA", h = d$d$h_CN, v = d$v, 
                         isKeptCnv = isKeptCnv, isMatchedCnv = isMatchedCnv, yaxt = "n", 
                         col = nrCol, cnvs = cnvs, chromI = zc_i)
        # }

        # plot CN vs. window index, i.e., post-normalization (color as final CN call)
        if(!short) plotCellByWindow(d$d$cww$CN, "CN corr.", h = d$d$h_CN, v = d$v, 
                         isKeptCnv = isKeptCnv, isMatchedCnv = isMatchedCnv, yaxt = "n", 
                         col = d$d$colNA_CN, cnvs = cnvs, chromI = zc_i, 
                         sampleHMM = d$d$sampleHMM, maskedHMM = d$d$maskedHMM)
    }
}
plotFAR <- function(cell, d, zoomChrom){
    par(mar = c(0.1, 0, 0.1, 0), cex = 1) 
    zc_i <- d$w[chrom == zoomChrom, i]
    y <- d$d$cww$FAR[zc_i]
    x <- 1:length(y)    
    plot(NA, NA, bty = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
         xlab = NULL, ylab = "FAR", xlim = range(x), ylim = c(-0.5, 2))
    abline(h = c(0, 0.5, 1), v = d$v, col = "grey")
    points(x, y, pch = 19, cex = 0.4, col = adjustColorAlpha(d$d$colNA_NR[zc_i], 0.75) )
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
getCellPlotPanel <- function(sourceId, project, cell, settings, level, clickAction,
                             widthChunks, layoutMatrix, zoomChrom = NULL, 
                             isKeptCnv = NULL, isMatchedCnv = NULL, force = FALSE, 
                             cnvs = NULL, short = FALSE, repPlot = FALSE){
    cellKey <- paste(sourceId, cell$cell_id)       
    d <- cellCompositeData[[cellKey]]
    pngFile <- getCellPngFile(project, cell, level, d)
    mustCreate <- force || !file.exists(pngFile)
    if(mustCreate) {
        png(
            filename = pngFile,
            width  = 1.5 * widthChunks, 
            height = 1.35 / (short + repPlot + 1), 
            units = "in", 
            pointsize = 7,
            bg = "white",
            res = 96, # i.e., optimized for on-screen display
            type = "cairo"
        )
        layout(layoutMatrix)
        if(repPlot) plotFAR(cell, d, zoomChrom)
        else plotCellQC(cell, d, zoomChrom, isKeptCnv, isMatchedCnv, cnvs, short)
        dev.off()
    }
    tags$img(
        src = pngFileToBase64(pngFile), 
        style = "vertical-align: top; cursor: pointer;", 
        class = "cellWindowsPlot",
        "data-source_id" = sourceId,
        "data-cell_id" = cell$cell_id,
        "data-action" = clickAction
    )
}
getCellCompositePlot <- function(sourceId, project, cell, settings){
    widthChunks <- 6
    force <- FALSE
    getCellPlotPanel(
        sourceId, project, cell, settings, 
        level = paste0("genome", ".tall"),         
        clickAction = "zoomChrom",
        widthChunks = genomeWidthChunks,
        layoutMatrix = matrix(c(c(1,1), rep(c(2,3), genomeWidthChunks - 1)), nrow = 2, ncol = genomeWidthChunks),
        force = force
    )  
}
getZoomChromPlot <- function(sourceId, project, cell, settings, zoomChrom, 
                             isKeptCnv = NULL, force = FALSE){
    if(is.null(zoomChrom)) return("")
    getCellPlotPanel( # two-row display of a single-cell on a chromosome, for use mainly by Keep/Reject
        sourceId, project, cell, settings, 
        level = paste0(zoomChrom, ".tall"),         
        clickAction = "keepReject", 
        widthChunks = zoomChromWidthChunks,
        layoutMatrix = matrix(rep(c(1,2), zoomChromWidthChunks), nrow = 2, ncol = zoomChromWidthChunks),
        zoomChrom = zoomChrom,
        isKeptCnv = isKeptCnv,
        force = force
    )  
}
getShortChromPlot <- function(sourceId, project, cell, settings, zoomChrom, cnvs, 
                              isMatchedCnv = NULL, force = FALSE){
    getCellPlotPanel(
        sourceId, project, cell, settings, 
        level = paste0(zoomChrom, ".short"),         
        clickAction = "toggleMatch", 
        widthChunks = matchChromWidthChunks,
        layoutMatrix = matrix(rep(1, matchChromWidthChunks), nrow = 1, ncol = matchChromWidthChunks),
        zoomChrom = zoomChrom,
        isMatchedCnv = isMatchedCnv,
        force = force,
        cnvs = cnvs,
        short = TRUE
    )  
}
getRepChromPlot <- function(sourceId, project, cell, settings, zoomChrom, force = FALSE){
    getCellPlotPanel(
        sourceId, project, cell, settings, 
        level = paste0(zoomChrom, ".replication"),         
        clickAction = "repClick", 
        widthChunks = matchChromWidthChunks,
        layoutMatrix = matrix(rep(1, matchChromWidthChunks), nrow = 1, ncol = matchChromWidthChunks),
        zoomChrom = zoomChrom,
        force = force,
        repPlot = TRUE
    ) 
}
setCellCompositeData <- function(sourceId, project, cell, settings, level,
                                 getReplicating = NULL, getKeep = NULL, force = FALSE){
    cellKey <- paste(sourceId, cell$cell_id)                            
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
    keyIsLoaded <- !is.null(cellCompositeData[[cellKey]]) &&
                   cellCompositeData[[cellKey]]$key == key
    dataIsLoaded <- keyIsLoaded && !is.null(cellCompositeData[[cellKey]]$d)
    if(!force && imageFileExists){
        if(!keyIsLoaded) cellCompositeData[[cellKey]] <<- list(
            key = key,
            shapeModel = shapeModel,
            replicationModel = replicationModel,
            col = divCol
        )
        return(cellKey)
    } else if(dataIsLoaded){
        return(cellKey)
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

        RPA_post_d <- cw_post$NR_wms[I_d] / (if(replicationModel$isSequential) cww$HMM[I_d] else cww$NA_[I_d])
        gc_fit <- if(replicationModel$isSequential) cww$gc_fit else cell$replicationModel[[shapeModel$key]]$gc_fit
        RPA_pre_w <- predict(gc_fit, w$gc_fraction, type = 'mu') # adjustedPeak
        cn_pre_w <- cw_pre$NR_wms / RPA_pre_w

        minWindows <- settings$get("Page_Options", "Minimum_Windows_Per_CNV", 20)
        cnvsType <- getCnvsType(settings)
        cnvs <- project$cnvs[[shapeModel$key]][type == cnvsType & cell_id == cell$cell_id & nWindows >= minWindows]
        cnvs <- lapply(seq_len(nrow(cnvs)), function(j){
            x <- cnvs[j]
            w_i <- w[chrom == x$chrom & start ==x$start, i]
            list(
                key = getCnvKey(x, sourceId),
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
            RPA_post_d = RPA_post_d,
            gc_fit = gc_fit,
            cn_pre_w = cn_pre_w,
            cnvs = cnvs,
            referenceCN = referenceCN,
            ploidy = cell$ploidy,
            sampleHMM = sampleHMM,
            maskedHMM = maskedHMM
        )
    }
    cellCompositeData[[cellKey]] <<- list(
        key = key,
        shapeModel = shapeModel,
        replicationModel = replicationModel,
        col = divCol,
        w = w,
        v = v,
        d = d
    )
    cellKey
}

#----------------------------------------------------------------------
# create a div of metadata about a cell and buttons to manipulate it
#----------------------------------------------------------------------
getCellSummary <- function(project, cell, buttons, short = FALSE){
    colData <- project$colData[cell$cell_id]
    desc <- project$manifest[as.character(Sample_ID) == cell$cell_id, Description]
    tags$div(
        style = if(short) "height: 64px;" else "",
        class = "cellSummaryPanel",
        tags$div(tags$strong( paste("cell:", colData$cell_id, '(', desc, ')') )), 
        tags$div(paste("windowPower: ", colData$windowPower, '(', commify(2 ** colData$windowPower * 20), 'kb )')),
        if(cell$badCell) tags$div("bad cell, not processed") else if(short) project$metadata$sample_id else tagList(
            tags$div(paste("repGcRatio:", round(colData$repGcRatio, 2), '(', colData$modelType, ',', round(colData$fractionS, 2), ')')), 
            tags$div(paste("cnsd:", round(colData$cnsd, 2), 
                           '(', if(colData$keep) "" else "not", "passed", ')',
                           '(', colData$sex, ', expect', colData$expectedSex, ')')),
            if(!is.null(buttons)) buttons(colData, cell) else "" 
        ) 
    )
}

#----------------------------------------------------------------------
# master function for plotting a single cell
#----------------------------------------------------------------------
cellCompositeData <- list()
createGenomeLabelRow <- function(project, cell){ # NOTE: don't make this dynamic to zoom, it forces an update of the entire genome panel
    chroms <- project$windows[[cell$windowPower + 1]][, .N, by = "chrom"]   
    chroms[, ":="(
        chrom = gsub("chr", "", chrom),
        pixelWidth = genomeXWidth * N / sum(N)
    )]
    tags$div(
        style = paste("border: 1px solid grey;"),
        class = "cellPlotsLabelWrapper",
        tagList(
            lapply(1:nrow(chroms), function(i) {
                chrom <- chroms[i, chrom]
                W <- chroms[i, pixelWidth]
                tags$div(
                    tags$strong(chrom), 
                    class = "chromLabelDiv", 
                    style = paste0("width: ", W, "px;", " max-width: ", W, "px;")
                )
            })
        )
    )
}
createZoomLabelRow <- function(session, zoomChrom, close = FALSE, commit = FALSE, leftPad = FALSE){
    tags$div(
        style = paste("border: 1px solid grey;", if(leftPad) "padding-left: 270px;" else ""),
        class = "cellPlotsLabelWrapper zoomPlotsLabelWrapper",
        tagList(
            tags$div(
                actionLink(session$ns("prevZoomChrom"), "previous", style = "float: left; margin-left: 10px;"),
                tags$strong(zoomChrom),                  
                if(close)  actionLink(session$ns("closeZoomPanel"), "close",  style = "margin-left: 10px;") else "",
                if(commit) tags$strong( actionLink(session$ns("commitCnvGroup"), "COMMIT GROUP", style = "margin-left: 10px;") ) else "",  
                actionLink(session$ns("nextZoomChrom"), "next", style = "float: right; margin-right: 10px;"),
                class = "chromLabelDiv zoomLabelDiv", 
                style = paste0("width: 100%;")
            )
        )
    )
}
plotOneCellUI_genome <- function(sourceId, project, cell, settings, buttons = NULL, 
                                 getReplicating = NULL, getKeep = NULL){
    req(project)
    force <- FALSE
    cellKey <- setCellCompositeData(sourceId, project, cell, settings, "genome", getReplicating, getKeep, force = force)
    tags$div(
        style = paste("border: 1px solid", cellCompositeData[[cellKey]]$col, ";"),
        class = "cellPlotWrapper",
        tagList(
            getCellSummary(project, cell, buttons),          
            getCellCompositePlot(sourceId, project, cell, settings)
        )
    )
}
plotOneCellUI_chrom <- function(sourceId, project, cell, settings, zoomChrom, invalidateZoomedCell,
                                 getReplicating = NULL, getKeep = NULL, isKeptCnv = NULL){
    req(project)
    force <- !is.null(invalidateZoomedCell) && cell$cell_id == invalidateZoomedCell
    cellKey <- setCellCompositeData(sourceId, project, cell, settings, zoomChrom, 
                         getReplicating, getKeep, force = force)
    tags$div(
        style = paste("border: 1px solid", cellCompositeData[[cellKey]]$col, ";"),
        class = "cellPlotWrapper",
        getZoomChromPlot(sourceId, project, cell, settings, zoomChrom, isKeptCnv = isKeptCnv, force = force)
    )
}
plotOneCellUI_match <- function(sourceId, project, cell, settings, zoomChrom, invalidateZoomedCells, 
                                isMatchedCnv, cnvs){
    req(project)
    force <- paste(sourceId, cell$cell_id) %in% invalidateZoomedCells

    force <- TRUE

    cellKey <- setCellCompositeData(sourceId, project, cell, settings, zoomChrom, force = force)
    tags$div(
        style = paste("height: 66px; border: 1px solid", cellCompositeData[[cellKey]]$col, ";"),
        class = "cellPlotWrapper",
        tagList(
            getCellSummary(project, cell, buttons = NULL, short = TRUE),
            getShortChromPlot(sourceId, project, cell, settings, zoomChrom, cnvs, 
                             isMatchedCnv = isMatchedCnv, force = force)
        )
    ) 
}
plotOneCellUI_FAR <- function(sourceId, project, cell, settings, zoomChrom){
    req(project)
    force = TRUE
    cellKey <- setCellCompositeData(sourceId, project, cell, settings, zoomChrom, force = force)
    tags$div(
        style = paste("height: 66px; border: 1px solid;"),
        class = "cellPlotWrapper",
        tagList(
            getRepChromPlot(sourceId, project, cell, settings, zoomChrom, force = force)
        )
    ) 
}

#----------------------------------------------------------------------
# handle plot click
#----------------------------------------------------------------------
genomeXOffset <- (1.5 + 0.5) * 96
genomeWidthChunks <- 7
genomeXWidth <- (1.5 * (genomeWidthChunks - 1) - 0.5) * 96
zoomChromXOffset <- 0
zoomChromWidthChunks <- 3
zoomChromXWidth <- 1.5 * zoomChromWidthChunks * 96
matchChromXOffset <- 0
matchChromWidthChunks <- 5
matchChromXWidth <- 1.5 * matchChromWidthChunks * 96
handleCellPlotClick <- function(click, zoomChrom, zoomTargetWindow){
    expandClick <- function(offset, width){
        x <- click$coord$x - offset
        project <- getScCnvProjectData(click$data$source_id)
        cell_id <- as.character(click$data$cell_id) # gets sent to us as an integer by Shiny.setInputValue
        cell <- project$cells[[cell_id]]
        w <- project$windows[[cell$windowPower + 1]] 
        list(
            x = x,
            w = w,
            cell_id = cell_id,
            project = project
        )       
    }
    getTargetWindow <- function(offset, width){
        click <- expandClick(offset, width)
        w_c <- click$w[chrom == zoomChrom()]
        n_w <- nrow(w_c)
        I <- round(click$x / width * n_w, 0)
        I <- pmin(n_w, pmax(1, I))
        w_c[I, ]
    }
    switch(
        click$data$action,
        zoomChrom = {
            click <- expandClick(genomeXOffset, genomeXWidth)
            chrom <- click$w[round(click$x / genomeXWidth * nrow(click$w), 0), chrom]
            zoomChrom(chrom)
        },
        keepReject = {
            zoomTargetWindow(list(
                targetWindow = getTargetWindow(zoomChromXOffset, zoomChromXWidth), 
                cell_id = as.character(click$data$cell_id),
                entropy = sample(1e8, 1) # ensure click happens even if mouse doesn't move
            ))
        },
        toggleMatch = {
            zoomTargetWindow(list(
                targetWindow = getTargetWindow(matchChromXOffset, matchChromXWidth), 
                click = click,
                entropy = sample(1e8, 1)
            ))
        }
    )
}
