# create additional informative QC plots
renderPlot <- function(plotFn, ..., suffix = NULL){
    filename <- if(is.null(suffix)) paste(plotFn, "png", sep = ".")
                else paste(plotFn, suffix, "png", sep = ".")
    pngFile <- file.path(env$PLOTS_DIR, filename)
    png(pngFile, width = 6, height = 6, units = "in", pointsize = 7, res = 600, type = "cairo")
    get(plotFn)(..., suffix)
    dev.off()
}
insertSize_vs_eventSize <- function(nodes, suffix, ...){
    eS <- nodes[, log10(ifelse(eventSize == 0, 250 * 1e6, eventSize))]
    iS <- nodes[, sapply(insertSize, function(x) if(x>0) log10(x) else if(x<0) -log10(-x) else 0)]
    plot(
        NA, 
        NA, 
        main = paste(suffix, "High-Quality Junctions"), 
        xlab = "Event Size (log10(bp))",
        ylab = "Insert Size (log10(bp))",
        xlim = range(eS, na.rm = TRUE),
        ylim = range(iS, na.rm = TRUE),
        col = NA
    )
    abline(
        v = 1:8,
        h = c(0, log10(c(1:10,100,1000)), -log10(c(1:10,100,1000))),
        col = rgb(0,0,0,0.1)
    )
    points(
        jitter(eS, amount = 0.1), 
        jitter(iS, amount = 0.1), 
        pch = 19, 
        cex = 0.25, 
        col = unlist(edgeTypeColors[nodes$edgeType])
    )
}
eventSizeDistribution <- function(nodes, suffix, ...){
    eS <- nodes[, log10(ifelse(eventSize == 0, 250 * 1e6, eventSize))]
    plot(
        density(eS, n = 512 * 4),
        main = paste(suffix, "High-Quality Junctions"), 
        xlab = "log10(eventSize)",
        ylab = "Density",
        xlim = log10(c(100, 250 * 1e6))
    )
}
fractionChimeric_vs_eventSize <- function(nodes, suffix, ...){
    eS <- nodes[, log10(ifelse(eventSize == 0, 250 * 1e6, eventSize))]
    plot(
        eS, 
        jitter(as.integer(nodes$isChimeric), amount = 0.1), 
        main = paste(suffix, "High-Quality Junctions"), 
        xlab = "Event Size (log10(bp))",
        ylab = "Fraction with an Adapter",
        xlim = log10(c(100, 250 * 1e6)),
        ylim = c(-0.1, 1.1),
        pch = 19, 
        cex = 0.25, 
        col = unlist(edgeTypeColors[nodes$edgeType])
    )
}
fractionChimeric_vs_insertSize <- function(nodes, suffix, ...){
    plot(
        jitter(nodes$insertSize,       amount = 0.5), 
        jitter(as.integer(nodes$isChimeric), amount = 0.1), 
        main = paste(suffix, "High-Quality Junctions"), 
        xlab = "Insertion Size (bp)",
        ylab = "Fraction with an Adapter",
        xlim = c(-100, 100),
        ylim = c(-0.1, 1.1),
        pch = 19, 
        cex = 0.25, 
        col = ifelse(nodes$edgeType == edgeTypes$TRANSLOCATION, rgb(0.8, 0, 0, 0.1), rgb(0, 0, 1, 0.1))
    )
}

fractionChimeric_vs_nInstances <- function(junctions, suffix, ...){
    plot(
        jitter(junctions$nInstances,       amount = 0.5), 
        jitter(junctions$fractionChimeric, amount = 0.1), 
        main = paste(suffix, "Unique High-Quality Junctions"),
        xlab = "# of Instances",
        ylab = "Fraction with an Adapter",
        xlim = c(0, 30),
        ylim = c(-0.1, 1.1),
        pch = 19, 
        cex = 0.25, 
        col = ifelse(junctions$edgeType == edgeTypes$TRANSLOCATION, rgb(0.8, 0, 0, 0.1), rgb(0, 0, 1, 0.1))
    )
}

adapterScore_vs_position <- function(trainingSet, nodes, endN, suffix){
    xcol <- if(endN == 5) "end" else "start"
    trainCol    <- paste0("trainable",  endN)
    xcol        <- paste0(xcol,         endN)
    ycol        <- paste0("score",      endN)
    predictCol  <- paste0("hasAdapter", endN)
    trainingSet <- trainingSet[trainingSet[[trainCol]] == TRUE]
    nodes <- nodes[sample(.N, min(1000, .N))]
    plot(
        jitter(trainingSet[[xcol]], amount = 0.5), 
        jitter(trainingSet[[ycol]], amount = 0.5), 
        main = if(is.null(suffix)) "Training Set plus 1K High-Quality SV Junctions" else suffix,
        xlab = xcol,
        ylab = ycol,
        pch = 19, cex = 0.25, 
        col = "grey"
    )
    points(
        jitter(nodes[[xcol]], amount = 0.5), 
        jitter(nodes[[ycol]], amount = 0.5), 
        pch = 19, cex = 0.25, 
        col = ifelse(nodes[[predictCol]] == TRUE, "red3", "blue")
    )
}



