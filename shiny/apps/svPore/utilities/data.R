# support for svPore data loading
svNodesReactive <- function(sourceId, session) reactive({ 
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading source nodes")
    x <- readRDS(getSourceFilePath(sourceId, "svNodesFile")) 
    stopSpinner(session)
    x[, isChrM := chrom1 == "chrM" | chrom2 == "chrM"]

    # x[1:100000]
    x

})
filteredNodesReactive <- function(svNodes, input) reactive({
    d <- svNodes()
    req(d)
    startSpinner(session, message = "filtering source nodes")        
    d <- d[edgeType != edgeTypes$ALIGNMENT & mapQ >= input$minMapQ]
    if(input$edgeType != "X") d <- d[edgeType == input$edgeType]
    if(!is.null(input$filters)){
        invert <- "invert" %in% input$filters
        if("passedBandwidth" %in% input$filters) d <- d[passedBandwidth == !invert]
        if("hasAdapter" %in% input$filters) d <- d[hasAdapter == invert]
        if("isChrM" %in% input$filters) d <- d[isChrM == invert]            
    }
    stopSpinner(session)
    req(nrow(d) > 0)
    d
})    
getSummaryValue <- function(d, column){
    switch(
        column,
        eventSize  = d[, log10(ifelse(eventSize == 0, 250 * 1e6, eventSize))],
        insertSize = d[, sapply(insertSize, function(x) if(x>0) log10(x) else if(x<0) -log10(-x) else 0)],
        d[[column]]
    )
}
getSummaryColor <- function(d, input){
    x <- d[[input$colorColumn]]
    edgeTypeColors <- list(
        "A" = rgb(0, 0, 0, input$opacity), # the single type for a contiguous aligned segment
        "T" = rgb(1, 0, 0, input$opacity), # edge/junction types (might be several per source molecule)
        "V" = rgb(0, 0.8, 0, input$opacity),
        "U" = rgb(0.8, 0.8, 0, input$opacity),
        "D" = rgb(0, 0, 1, input$opacity),
        "?" = rgb(0, 0, 0, input$opacity),
        "I" = rgb(0, 0, 0, input$opacity), 
        "M" = rgb(0, 0, 0, input$opacity),
        "P" = rgb(0, 0, 0, input$opacity),
        "F" = rgb(0, 0, 0, input$opacity),
        "R" = rgb(0, 0, 0, input$opacity),
        "Q" = rgb(0, 0, 0, input$opacity)
    )
    colorLogical <- function(xx, opacity) ifelse(xx, rgb(0.8, 0, 0, opacity), rgb(0, 0, 0, opacity))
    switch(
        input$colorColumn,
        edgeType = unlist(edgeTypeColors[x]),
        passedBandwidth = colorLogical(!x, input$opacity),
        hasAdapter = colorLogical(x, input$opacity),
        isChrM = colorLogical(x, input$opacity),
        rgb(0, 0, 0, input$opacity)
    )
}
plotNodesReactive <- function(filteredNodes, input) reactive({
    d <- filteredNodes()
    req(d)
    startSpinner(session, message = "parsing plot data")
    d <- data.table(
        qName = d$qName,
        junction = d$junction,
        x = getSummaryValue(d, input$xAxisColumn),
        y = getSummaryValue(d, input$yAxisColumn),
        color = getSummaryColor(d, input)
    )
    d[, ":="(
        xj = jitter(as.numeric(x), amount = 0.1),
        yj = jitter(as.numeric(y), amount = 0.1)
    )]
    stopSpinner(session)
    d
})    
#  [1] "junction"         "qName"            "edge"             "node1" "cigar1"
#  [5] "node2" "cigar2"            "edgeType"         "mapQ"             "eventSize"       
#  [9] "insertSize"       "xStart"           "xEnd"             "edgeClass"       
# [13] "nStrands"         "chromIndex1"      "chrom1"           "windowIndex1"    
# [17] "strand1"          "chromIndex2"      "chrom2"           "windowIndex2"    
# [21] "strand2"          "nInstances"       "nMolecules"       "passedBandwidth" 
# [25] "score3"           "score5"           "start3"           "end5"
# [29] "hasAdapter3"      "hasAdapter5"      "hasAdapter"       "fractionChimeric"
# [33] "segment"          "segmentName"
segmentsReactive <- function(svNodes, plotNodes, summaryPlot) reactive({
    svNodes <- svNodes()
    plotNodes <- plotNodes()
    coord <- summaryPlot$click()$coord
    distance <- plotNodes[, sqrt((x - coord$x)**2 + (y - coord$y)**2)]
    indexNode <- plotNodes[which.min(distance)]
    indexJunction <- indexNode[, junction]
    indexQName <- indexNode[, qName]
    qNames <- svNodes[junction == indexJunction & qName != indexQName, unique(qName)]
    if(length(qNames) > 2) qNames <- sample(qNames, 2)
    qNames <- c(indexQName, qNames)
    x <- svNodes[qName %in% qNames][order(qName, edge)]
    
    padding <- 10000
    maxI <- nrow(x)
    maxQSize <- x[edgeType == edgeTypes$ALIGNMENT, max(xEnd)]
    x[, locus := as.integer(NA)]
    getAlnRPos <- function(qName_, edge_) {
        xq <- x[qName == qName_]
        if(edge_ == 1) xq[edge_ + 1, xStart] else xq[edge_ - 1, xEnd] 
    }
    nLoci <- 0
    fillLocus <- function(x){
        pending <- x[, edgeType == edgeTypes$ALIGNMENT & is.na(locus)]
        i <- x[, min(which(pending), na.rm = TRUE)]
        chrom <- x[i, chrom1]
        rPos  <- x[i, getAlnRPos(qName, edge)]
        rRange <- c(
            rPos - padding - maxQSize,
            rPos + padding + maxQSize 
        )
        x[, locus := {
            if(!is.na(locus) || edgeType != edgeTypes$ALIGNMENT) locus
            else if(chrom1 == chrom && getAlnRPos(qName, edge) %between% rRange) abs(x[i, node1])
            else as.integer(NA)   
        }, by = .(qName, edge)]
        nLoci <<- nLoci + 1
        x
    }
    while(x[edgeType == edgeTypes$ALIGNMENT, any(is.na(locus))]) x <- fillLocus(x)

    dt <- do.call(rbind, lapply(1:maxI, function(i){
        if(x[i, edgeType == edgeTypes$ALIGNMENT]) data.table(
            qName = x[i, qName],
            edge  = x[i, edge],
            edgeType = edgeTypes$ALIGNMENT,
            cigar   = x[i, cigar1],
            locus1  = x[i, locus],
            locus2  = x[i, locus],
            chrom   = x[i, chrom1],
            strand  = x[i, sign(node1)],
            qStart  = x[i, xStart],
            qEnd    = x[i, xEnd],
            rStart  = if(x[i, edge > 1])    x[i - 1, xEnd]   else getAlnStartPos(x[i + 1, xStart], x[i, node1], x[i, cigar1]),
            rEnd    = if(x[i, edge] < x[, max(edge)]) x[i + 1, xStart] else   getAlnEndPos(x[i - 1, xEnd],   x[i, node1], x[i, cigar1]),
            qOffset = as.integer(NA),
            rOffset = as.integer(NA),
            label   = x[i, paste0("MAPQ: ", mapQ, "; bp:", commify(eventSize))]
        ) else data.table(
            qName = x[i, qName],
            edge  = x[i, edge],
            edgeType = x[i, edgeType],
            cigar    = as.character(NA),
            locus1  = x[i - 1, locus],
            locus2  = as.integer(NA),
            chrom   = as.character(NA),
            strand  = as.integer(NA),
            qStart  = x[i - 1, xEnd],
            qEnd    = x[i + 1, xStart],
            rStart  = x[i, xStart],
            rEnd    = x[i, xEnd],
            qOffset = x[i, insertSize],
            rOffset = x[i, eventSize],
            label   = as.character(NA)
        )
    }))
    for(i in 1:maxI) if(dt[i, edgeType != edgeTypes$ALIGNMENT]){
        ln2 <- dt[i + 1, locus1]
        dt[i, locus2 := ln2]
    }
    list(
        qNames = qNames,
        loci = dt[edgeType == edgeTypes$ALIGNMENT, locus1],
        uniqueLoci = dt[edgeType == edgeTypes$ALIGNMENT, unique(locus1)],
        nLoci = nLoci,
        dt = dt,
        maxI = maxI
    )
})
