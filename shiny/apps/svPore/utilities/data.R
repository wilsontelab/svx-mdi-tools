# support for svPore data loading
svEdgesReactive <- function(sourceId, session) reactive({ 
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading source edges")
    x <- readRDS(getSourceFilePath(sourceId, "svEdgesFile")) 
    stopSpinner(session)
    x[, isChrM := chrom1 == "chrM" | chrom2 == "chrM"]
    # x[1:100000]
    setkey(x, qName, blockN, edgeN)
    x
})
# Classes ‘data.table’ and 'data.frame':  56116 obs. of  57 variables:
#  $ qName                : chr  "00064a48-c298-4bb8-a1de-4fc12ddc9510" "00064a48-
# c298-4bb8-a1de-4fc12ddc9510" "00064a48-c298-4bb8-a1de-4fc12ddc9510" "0008dd97-96
# 67-4819-aed9-e1d31f68ba42" ...
#  $ junctionKey          : chr  NA "13:1:1:1:92730978:29844058:0" NA NA ...
#  $ blockN               : int  1 2 3 1 2 3 1 2 3 1 ...
#  $ edgeN                : int  1 2 3 1 2 3 1 2 3 1 ...
#  $ node1                : int  -29867 -29844 -2193327 290649 290667 -1545388 -98
# 0018 -980018 2506843 -1285732 ...
#  $ blastIdentity        : num  0.99 0.99 0.994 0.991 0.991 ...
#  $ gapCompressedIdentity: num  0.993 0.993 0.994 0.994 0.994 ...
#  $ node2                : int  -29844 -2193327 -2193324 290667 -1545388 -1545370
#  -980018 2506843 2506843 -1285731 ...
#  $ edgeType             : chr  "A" "T" "A" "A" ...
#  $ mapQ                 : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ eventSize            : int  23245 0 2986 17920 0 18304 305 0 492 1654 ...    
#  $ insertSize           : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ xStart               : int  29 29844058 23260 38 42279419 17916 35 94252532 3
# 41 26 ...
#  $ xEnd                 : int  23260 92730978 26246 17916 144881594 36087 341 91
# 763031 836 1656 ...
#  $ edgeClass            : int  NA 2 NA NA 2 NA NA 2 NA NA ...
#  $ nStrands             : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ chrom1               : chr  "chr1" "chr1" "chr13" "chr2" ...
#  $ chromIndex1          : int  1 1 13 2 2 8 5 5 16 7 ...
#  $ windowIndex1         : num  29867 29844 92730 42261 42279 ...
#  $ genomeIndex1         : int  29867 29844 2193327 290649 290667 1545388 980018 
# 980018 2506843 1285732 ...
#  $ strand1              : chr  "-" "-" "-" "+" ...
#  $ chrom2               : chr  "chr1" "chr13" "chr13" "chr2" ...
#  $ chromIndex2          : int  1 13 13 2 8 8 5 16 16 7 ...
#  $ windowIndex2         : num  29844 92730 92727 42279 144881 ...
#  $ genomeIndex2         : int  29844 2193327 2193324 290667 1545388 1545370 9800
# 18 2506843 2506843 1285731 ...
#  $ strand2              : chr  "-" "-" "-" "+" ...
#  $ edgeId               : int  12 13 14 21 22 23 24 25 26 27 ...
#  $ nEdges               : int  3 3 3 3 3 3 3 3 3 3 ...
#  $ passedFlankCheck     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ passedBandwidth      : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ score3               : num  NA 4.5 NA NA 4 NA NA 5 NA NA ...
#  $ score5               : num  NA 4.5 NA NA 4 NA NA 5 NA NA ...
#  $ start3               : num  NA 18 NA NA -2 NA NA 22 NA NA ...
#  $ end5                 : num  NA 4 NA NA -21 NA NA -27 NA NA ...
#  $ hasAdapter3          : logi  NA FALSE NA NA FALSE NA ...
#  $ hasAdapter5          : logi  NA FALSE NA NA FALSE NA ...
#  $ hasAdapter           : logi  NA FALSE NA NA FALSE NA ...
#  $ isCanonical          : logi  NA FALSE NA NA TRUE NA ...
#  $ cChromIndex1         : int  1 13 13 2 2 8 5 16 16 7 ...
#  $ cChromIndex2         : int  1 1 13 2 8 8 5 5 16 7 ...
#  $ cStrand1             : int  NA 1 NA NA 1 NA NA -1 NA NA ...
#  $ cStrand2             : int  NA 1 NA NA -1 NA NA 1 NA NA ...
#  $ cRefPos1             : int  NA 92730978 NA NA 42279419 NA NA 91763031 NA NA .
# ..
#  $ cRefPos2             : int  NA 29844058 NA NA 144881594 NA NA 94252532 NA NA 
# ...
#  $ indexJunctionKey     : chr  NA "13:1:1:1:92730978:29844058:0" NA NA ...      
#  $ isIndexJunction      : logi  NA TRUE NA NA TRUE NA ...
#  $ nCanonical           : int  NA 0 NA NA 1 NA NA 0 NA NA ...
#  $ nNonCanonical        : int  NA 1 NA NA 0 NA NA 1 NA NA ...
#  $ icRefPos1            : int  NA 92730978 NA NA 42279419 NA NA 91763031 NA NA .
# ..
#  $ icRefPos2            : int  NA 29844058 NA NA 144881594 NA NA 94252532 NA NA 
# ...
#  $ iInsertSize          : int  NA 0 NA NA 0 NA NA 0 NA NA ...
#  $ iRefPos1             : int  NA 29844058 NA NA 42279419 NA NA 94252532 NA NA .
# ..
#  $ iRefPos2             : int  NA 92730978 NA NA 144881594 NA NA 91763031 NA NA 
# ...
#  $ icPathKey            : chr  "2193324:::2193322474:0:29844057:::29867" "219332
# 4:::2193322474:0:29844057:::29867" "2193324:::2193322474:0:29844057:::29867" "29
# 0649:::290666746:0:-1545386061:::-1545370" ...
#  $ isDuplex             : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ retained             : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#  $ nJunctionInstances   : int  NA 1 NA NA 1 NA NA 1 NA NA ...
#  - attr(*, ".internal.selfref")=<externalptr>
#  - attr(*, "sorted")= chr [1:3] "qName" "blockN" "edgeN"
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
    svEdges <- svEdges()
    filteredJunctions <- filteredJunctions()
    plotJunctions <- plotJunctions()
    coord <- summaryPlot$click()$coord
    req(nrow(filteredJunctions) > 0, nrow(plotJunctions) > 0, coord, coord$x, coord$y)
    distance <- plotJunctions[, sqrt((x - coord$x)**2 + (y - coord$y)**2)]
    indexNode <- plotJunctions[which.min(distance)]
    indexJunctionKey_ <- indexNode[, if(is.na(indexJunctionKey)) "UNMATCHABLE_JUNCTION" else indexJunctionKey]
    indexQName <- indexNode[, qName]
    otherQNames <- filteredJunctions[indexJunctionKey == indexJunctionKey_ & qName != indexQName, unique(qName)]
    if(length(otherQNames) > 2) otherQNames <- sample(otherQNames, 2)
    qNames <- c(indexQName, otherQNames)
    x <- svEdges[qNames]
    padding <- 10000
    maxI <- nrow(x)
    maxQSize <- x[edgeType == edgeTypes$ALIGNMENT, max(xEnd)]
    x[, locus := NA_integer_]
    getAlnRPos <- function(qName_, edgeN_) {
        xq <- x[qName == qName_]
        if(edgeN_ == 1) xq[edgeN_ + 1, xStart] else xq[edgeN_ - 1, xEnd] 
    }
    nLoci <- 0
    fillLocus <- function(x){
        pending <- x[, edgeType == edgeTypes$ALIGNMENT & is.na(locus)]
        i <- x[, min(which(pending), na.rm = TRUE)]
        chrom <- x[i, chrom1]
        rPos  <- x[i, getAlnRPos(qName, edgeN)]
        rRange <- c(
            rPos - padding - maxQSize,
            rPos + padding + maxQSize 
        )
        x[, locus := {
            if(!is.na(locus) || edgeType != edgeTypes$ALIGNMENT) locus
            else if(chrom1 == chrom && getAlnRPos(qName, edgeN) %between% rRange) abs(x[i, node1])
            else as.integer(NA)   
        }, by = .(qName, edgeN)]
        nLoci <<- nLoci + 1
        x
    }
    while(x[edgeType == edgeTypes$ALIGNMENT, any(is.na(locus))]) x <- fillLocus(x)
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
