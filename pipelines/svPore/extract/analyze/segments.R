# segments are arbitrary length vectors of high quality nodes, including the outer molecule endpoints
# that describe a single, continuous, non-chimeric read

# concatenate nodes and edges into a single row per segment
collapseSegments <- function(edges){
    segments <- edges[, .(
        segmentLength = sum(eventSize[edgeType == edgeTypes$ALIGNMENT]),
        chroms = paste(unique(c(chrom1, chrom2)), collapse = ","),
        nJxns = sum(edgeType != edgeTypes$ALIGNMENT),
        nKeptJxns = sum(keptJunction),
        pathType = paste0(ifelse(edgeType == edgeTypes$ALIGNMENT | keptJunction, edgeType, tolower(edgeType)), collapse = ""),
        isClosedPath = chromIndex1[1] == chromIndex2[.N] && strand1[1] == strand2[.N],
        outerNode1 = rNode1[1], # original nodes of outer molecule endpoints, on strand(s) as sequenced
        nodePath = list(c(rbind(rNode1[keptJunction], rNode2[keptJunction]))), # reference junction node sequence, on strand(s) as sequenced
        outerNode2 = rNode2[.N],
        matchingSegmentsTmp = { # a matching segment has at least one fuzzy junction in common with the query segment
            x <- unique(unlist(matchingSegments))
            list(x[x != segmentName])
        }
    ), by = .(segmentName)]
    segments[, nMatchingSegments := sapply(matchingSegmentsTmp, function(x) length(unlist(x)))]
    segments
}
findMatchingSegments <- function(segments){
    I <- segments[, which(nKeptJxns > 0 & nMatchingSegments > 0)]
    if(length(I) == 0) return(data.table(segment1 = character(), segment2 = character(), N = integer()))
    x <- segments[I, .(
        matchingSegment = unlist(matchingSegmentsTmp)
    ), by = .(segmentName)]
    if(nrow(x) == 0) return(data.table(segment1 = character(), segment2 = character()))
    x[, {
        x <- sort(c(segmentName[1], matchingSegment[1])) # order the segment identifiers (not the nodes within them)
        .(
            segment1 = x[1],
            segment2 = x[2] 
        )
    }, by = .(segmentName, matchingSegment)][,
        .N,
        by = .(segment1, segment2) # two segments found to share at least one junction in common, thus are overlapping
    ]
}
analyzeSegmentsNetwork <- function(segments, segmentMatches){ # use igraph to identify and examine clusters of overlapping segments
    g <- graph_from_data_frame(segmentMatches, directed = FALSE)
    cmp <- components(g)
    I <- segments[, which(nKeptJxns > 0 & nMatchingSegments > 0)]
    clusters <- merge( # scan the network of segment matches for clusters of molecules sharing junctions
        data.table(segmentName = names(cmp$membership), cluster = cmp$membership),
        segments[I, segmentLength, by = .(segmentName)],
        by = "segmentName",
        all.x = TRUE
    )
    setkey(clusters, cluster, segmentLength)
    clusters[, { # assign a reference segment for each cluster as the the longest read in the group
        .(       # in the first round, that segment might well have an artifact junction that is inappropriately fusing two SV segments
            segmentName = segmentName,
            isReferenceSegment = segmentName == segmentName[.N],        
            refSegment = segmentName[.N],
            refSegmentLength = segmentLength[.N]
        )
    }, by = .(cluster)][,
        .SD,
        .SDcols = c("segmentName","isReferenceSegment","refSegment","refSegmentLength")
    ]
}

# there is always an even number of nodes as passed to these functions (2 per junction, 2 outer alignment nodes)
pathMatchTypes <- list(
    NONE        = 0,
    IDENTICAL   = 1, # IDENTICAL, DUPLEX, and SAME_NODES demand an equal number of nodes in the two paths
    DUPLEX      = 2, # IDENTICAL and DUPLEX require matching outer endpoints, i.e., molecules are considered to be the same source DNA molecule
    SAME_NODES  = 3, # SAME_NODES, NESTED, and OVERLAP do NOT consider molecule outer endpoints, i.e., molecules are considered to be different source DNA molecules
    NESTED      = 4, # NESTED has one segment entirely contained in another # TODO: not currently implemented
    OVERLAP     = 5  # OVERLAP has flanks on either side of each molecule
)
getOverlapMatch <- function(path1, path2){ # paths are arbitrary length vectors of high quality junction nodes, NOT including the outer molecule endpoints
    if(length(path2) > length(path1)){
        tmp <- path1
        path1 <- path2 # path1 is the longer of the two
        path2 <- tmp
    } 
    np2 <- length(path2)
    flank <- rep(NA_integer_, np2 - 2) 
    x1  <- c(flank,      path1,  flank)
    x1r <- c(flank, -rev(path1), flank)
    npos <- length(x1)
    for(register in 0:(npos - np2)){ 
        x2 <- c(
            rep(NA_integer_, register), 
            path2, 
            rep(NA_integer_, npos - np2 - register)
        )
        i <- !is.na(x1 + x2)
        if(x1[i] %is_exact_match% x2[i] || x1r[i] %is_exact_match% x2[i]) return(pathMatchTypes$OVERLAP)
    }
    pathMatchTypes$NONE
}
getSegmentMatch <- function(segment1, segment2){ # segment must be c(outerNode, junction, ..., outerNode)
    path1 <- segment1[2:(length(segment1) - 1)] # these junction nodes are reference nodes, support exact matching
    path2 <- segment2[2:(length(segment2) - 1)]
    if(length(segment1) == length(segment2)){
        if(path1 %is_exact_rev_match% path2 || path1 %is_exact_match% path2) {
            outer1 <- segment1[c(1, length(segment1))] # these are still initial nodes as aligned, require fuzzy matching
            outer2 <- segment2[c(1, length(segment2))]
            if(outer1 %is_fuzzy_rev_match% outer2) pathMatchTypes$DUPLEX 
            else if(outer1 %is_fuzzy_match% outer2) pathMatchTypes$IDENTICAL            
            else pathMatchTypes$SAME_NODES
        } else if(length(path1) > 2 || length(path2) > 2) getOverlapMatch(path1, path2)
        else pathMatchTypes$NONE
    } else getOverlapMatch(path1, path2)
}

# examine segment pairs for identity relationships, either on same or opposite strands
getOuterEndpointMatch <- function(outerNodes1, outerNodes2){ # these junction nodes are reference nodes, support exact matching
    if(outerNodes1 %is_fuzzy_rev_match% outerNodes2) return(pathMatchTypes$DUPLEX)
    if(outerNodes1 %is_fuzzy_match% outerNodes2)     return(pathMatchTypes$IDENTICAL)
    pathMatchTypes$NONE
}
getNodePathMatch <- function(nodePath1, nodePath2){ # these junction nodes are reference nodes, support exact matching
    if(length(nodePath1) != length(nodePath2))   return(pathMatchTypes$NONE)
    if(nodePath1 %is_exact_rev_match% nodePath2) return(pathMatchTypes$DUPLEX)
    if(nodePath1 %is_exact_match% nodePath2)     return(pathMatchTypes$IDENTICAL)
    pathMatchTypes$NONE
}
getSegmentMatchTypes <- function(segmentName1, nodePath1, outerNodes1, matchingSegmentNames){
    if(is.null(matchingSegmentNames)) return(
        data.table(segmentName1 = character(), segmentName2 = character(), 
                   nodePathMatchType = integer(), outerEndpointMatchType = integer())
    )
    segmentsToMatch[
        segmentName %in% matchingSegmentNames, # the list of segments identified above as having >=1 junction match to segment1
        .( 
            nodePathMatchType = getNodePathMatch(nodePath1, unlist(nodePath)),
            outerEndpointMatchType = getOuterEndpointMatch(outerNodes1, c(outerNode1, outerNode2))
        ), 
        by = .(segmentName)
    ][, .(
        segmentName1 = segmentName1,
        segmentName2 = segmentName,
        nodePathMatchType = nodePathMatchType,
        outerEndpointMatchType = outerEndpointMatchType
        # matchingSegments = list(paste(segmentName, matchType, sep = ":"))
    )]
}
scoreSegmentMatches <- function(segments){
    do.call(rbind, mclapply(1:nrow(segmentsToMatch), function(i){
        segmentsToMatch[i, getSegmentMatchTypes(segmentName, unlist(nodePath), c(outerNode1, outerNode2), unlist(matchingSegmentsTmp))]
    }, mc.cores = env$N_CPU))    
}
