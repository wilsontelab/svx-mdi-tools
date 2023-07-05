
# support for svPore data loading
junctionEdgesReactive <- function(sourceId, session) reactive({ 
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading junction edges")
    x <- readRDS(getSourceFilePath(sourceId, "edgesFile")) 
    stopSpinner(session)
    # x[, isChrM := chrom1 == "chrM" | chrom2 == "chrM"]
    # # x[1:100000]
    # setkey(x, qName, blockN, edgeN)
    x
})

junctionClustersReactive <- function(sourceId, session) reactive({ 
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading junction clusters")
    x <- readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
    stopSpinner(session)
    # x[, isChrM := chrom1 == "chrM" | chrom2 == "chrM"]
    # # x[1:100000]
    # setkey(x, qName, blockN, edgeN)
    data.table(a=1:10)    
})
# <<< stopSpinner: 
# Classes ‘data.table’ and 'data.frame':  86602 obs. of  22 variables:
#  $ clusterN             : num  1 2 3 4 5 6 7 8 9 10 ...
#  $ edgeType             : chr  "V" "V" "D" "D" ...
#  $ eventSize            : int  117668262 8692 98319533 74228034 6259 14735 0 0 0
#  0 ...
#  $ cChromIndex1         : int  10 10 10 10 10 10 10 10 10 10 ...
#  $ cChromIndex2         : int  10 10 10 10 10 10 11 12 12 12 ...
#  $ cStrand1             : int  -1 1 1 1 1 1 1 1 1 1 ...
#  $ cStrand2             : int  1 -1 1 1 1 1 -1 -1 -1 -1 ...
#  $ cRefPos1             : int  127004201 37137222 17203466 23997033 39328322 394
# 46704 1945052 116033555 78740899 97254217 ...
#  $ cRefPos2             : int  9335939 37145914 115522999 98225067 39334581 3943
# 1969 12149997 31364182 12102774 50310792 ...
#  $ insertSize           : int  -520 -2726 3 0 -4436 -4332 0 -6 -506 0 ...       
#  $ mapQ                 : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ gapCompressedIdentity: num  0.988 0.991 0.992 0.992 0.965 ...
#  $ baseQual             : num  41.7 36.6 NA NA 32.9 37.2 NA 50 40.9 NA ...      
#  $ alnBaseQual          : num  37 35.2 35.9 38.8 33.3 ...
#  $ alnSize              : int  30236 7979 3206 24925 17617 13893 2683 6026 3002 
# 1511 ...
#  $ samples              : chr  "GM24385_HG002" "GM24385_HG002" "GM24385_HG002" "
# GM24385_HG002" ...
#  $ nSamples             : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ nInstances           : int  2 2 2 2 2 4 2 2 2 2 ...
#  $ nCanonical           : int  1 2 1 1 0 1 1 1 1 1 ...
#  $ nNonCanonical        : int  1 0 1 1 2 3 1 1 1 1 ...
#  $ GM24385_HG002        : int  2 2 2 2 2 4 2 2 2 2 ...
#  $ segmentClusterN      : num  NA 106 NA NA 108 25 NA NA NA NA ...
#  - attr(*, ".internal.selfref")=<extern

filteredJunctionsReactive <- function(svEdges, input) reactive({
    d <- svEdges()
    req(d)
    startSpinner(session, message = "filtering source edges")  
    if(input$edgeType != "X") d <- d[edgeType == input$edgeType]
    if(input$showChrM == "No") d <- d[!(chrom1 == "chrM" | chrom2 == "chrM")]
    filter <- get(paste0("get", input$filterType))(d)
    d <- d[filter]
    message()
    print(d[, .N])
    print(d[, length(unique(indexJunctionKey))])
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
        passedFlankCheck = colorLogical(!x, input$opacity),
        passedBandwidth = colorLogical(!x, input$opacity),
        hasAdapter = colorLogical(x, input$opacity),
        isChrM = colorLogical(x, input$opacity),
        rgb(0, 0, 0, input$opacity)
    )
}
plotJunctionsReactive <- function(filteredJunctions, input) reactive({
    d <- filteredJunctions()
    req(d)
    startSpinner(session, message = "parsing plot data")
    d <- data.table(
        qName = d$qName,
        indexJunctionKey = d$indexJunctionKey,
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
segmentsReactive <- function(svEdges, filteredJunctions, plotJunctions, summaryPlot) reactive({
    # svEdges <- svEdges()
    # filteredJunctions <- filteredJunctions()
    # plotJunctions <- plotJunctions()
    # coord <- summaryPlot$click()$coord
    # req(nrow(filteredJunctions) > 0, nrow(plotJunctions) > 0, coord, coord$x, coord$y)
    # distance <- plotJunctions[, sqrt((x - coord$x)**2 + (y - coord$y)**2)]
    # indexNode <- plotJunctions[which.min(distance)]
    # indexJunctionKey_ <- indexNode[, if(is.na(indexJunctionKey)) "UNMATCHABLE_JUNCTION" else indexJunctionKey]
    # indexQName <- indexNode[, qName]
    # otherQNames <- filteredJunctions[indexJunctionKey == indexJunctionKey_ & qName != indexQName, unique(qName)]
    # if(length(otherQNames) > 2) otherQNames <- sample(otherQNames, 2)
    # qNames <- c(indexQName, otherQNames)
    # x <- svEdges[qNames]



    maxI <- nrow(x)
    dt <- do.call(rbind, lapply(1:maxI, function(i){
        if(x[i, edgeType == edgeTypes$ALIGNMENT]) data.table(
            qName = x[i, qName],
            edgeN  = x[i, edgeN],
            edgeType = edgeTypes$ALIGNMENT,
            cigar   = x[i, cigar],
            locus1  = x[i, locus],
            locus2  = x[i, locus],
            chrom   = x[i, chrom1],
            strand  = x[i, sign(node1)],
            qStart  = x[i, xStart],
            qEnd    = x[i, xEnd],
            rStart  = if(x[i, edgeN > 1])    x[i - 1, xEnd]   else getAlnStartPos(x[i + 1, xStart], x[i, node1], x[i, cigar]),
            rEnd    = if(x[i, edgeN] < x[, max(edgeN)]) x[i + 1, xStart] else   getAlnEndPos(x[i - 1, xEnd],   x[i, node1], x[i, cigar]),
            qOffset = as.integer(NA),
            rOffset = as.integer(NA),
            label   = x[i, paste0("MAPQ: ", mapQ, "; bp:", commify(eventSize))]
        ) else data.table(
            qName = x[i, qName],
            edgeN  = x[i, edgeN],
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

