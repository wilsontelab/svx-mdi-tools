#-------------------------------------------------------------------------------------
# segments are arbitrary length vectors of high quality nodes, including the outer molecule endpoints
# that describe a single, contiguous, non-chimeric source DNA span
#-------------------------------------------------------------------------------------

# split reads into potentially multiple segments, breaking reads at any untrustworthy/unvalidated junctions
getSplitSegmentN <- function(N, splitAt){
    junctionI <- seq(2, N, 2)
    if(!any(splitAt[junctionI])) 1 
    else if(N == 3) c(1, NA, 2)
    else {
        segmentNs <- 1
        for(i in junctionI){
            segmentNs <- c(
                segmentNs, 
                if(splitAt[i]) c(NA, segmentNs[i - 1] + 1)
                else rep(segmentNs[i - 1], 2)
            )
        }
        segmentNs
    }    
}
splitReadsToSegments <- function(edges){
    message("splitting reads at low-quality, chimeric, or unconfirmed junctions")
    fusable <- getFusableJunctions(edges)
    edges[, segmentN := getSplitSegmentN(.N, !fusable[.I]), by = .(qName)]
    edges <- edges[!is.na(segmentN)] # drop the false/untrusted junctions
    setkey(edges, qName, segmentN)

    # update blockN and edgeN within each split segment
    edges[, ":="(
        segmentName = paste(qName, segmentN, sep = "-"),   
        blockN = blockN - min(blockN) + 1,
        edgeN  = edgeN  - min(edgeN)  + 1
    ), by = .(qName, segmentN)]

    # remove split segments that do not contain any potential SV junctions
    # e.g., in A[T]A[T]A or A[T]AdA[T]A (where "d" was rejected but not split)
    # TODO: save these segments as alignments for haplotype assembly?
    matchable <- getMatchableJunctions(edges)
    x <- edges[, .(hasCandidateSvs = any(matchable[.I])), by = .(segmentName)]
    edges <- edges[segmentName %in% x[hasCandidateSvs == TRUE, segmentName]]
    setkey(edges, segmentName, blockN, edgeN)

    # prepare for assembly segment matching
    matchable <- getMatchableJunctions(edges)
    edges[matchable, ":="( 
        matchingSegments = list(unique(segmentName))
        # ,
        # nInstances = .N
    ), by = .(refJunctionKey)]
    edges
}

# concatenate nodes and edges into a single row per segment
collapseSegments <- function(edges){
    message("collapsing edges into contiguous, trusted read segments")
    isAlignment <- getAlignmentEdges(edges)
    isMatchable <- getMatchableJunctions(edges)
    segments <- edges[, {
        alignment <- isAlignment[.I]
        matchable <- isMatchable[.I]
        .(
            segmentLength = sum(eventSize[alignment]),
            chroms = paste(unique(c(chrom1, chrom2)), collapse = ","),
            nJxns = sum(!alignment),
            nMatchableJxns = sum(matchable),
            pathType = paste0(ifelse(alignment | matchable, edgeType, tolower(edgeType)), collapse = ""),
            isClosedPath = chromIndex1[1] == chromIndex2[.N] && strand1[1] == strand2[.N],
            outerNode1 = node1[1], # original nodes of outer molecule endpoints, on strand(s) as sequenced
            rNodePath = list(c(rbind(rNode1[matchable], rNode2[matchable]))), # reference junction node sequence, on strand(s) as sequenced
            outerNode2 = node2[.N],
            matchingSegmentsTmp = { # a matching segment has at least one fuzzy junction in common with the query segment
                x <- unique(unlist(matchingSegments))
                list(x[x != segmentName])
            }
        )
    }, by = .(segmentName)]
    segments[, ":="(
        nMatchingSegments = sapply(matchingSegmentsTmp, function(x) length(unlist(x))),
        i = 1:.N
    )]
    segments[, ":="( # strand orientation for segments is determined by the outermost alignments
        isCanonical = isCanonicalStrand(c(outerNode1, outerNode2))
    ), by = .(i)]
    segments[, ":="(
        cOuterNode1 = if(isCanonical) outerNode1 else -outerNode2,
        rcNodePath  = if(isCanonical) rNodePath  else list(-rev(unlist(rNodePath))),
        cOuterNode2 = if(isCanonical) outerNode2 else -outerNode1
    ), by = .(i)]
    segments[, ":="( # use exact matching on outer nodes for duplex purging
        pathKey = paste(cOuterNode1, unlist(rcNodePath), cOuterNode2, sep = ":", collapse = ":")
    ), by = .(i)]
    segments
}

# assemble the networks of segments with shared junctions, i.e., overlapping segments on a haplotype
findMatchingSegments <- function(segments){
    message("finding groups of overlapping segments by virtue of shared individual junctions")
    matchable <- getMatchableSegments(segments)
    nullTable <- data.table(segment1 = character(), segment2 = character(), N = integer())
    if(sum(matchable) == 0) return(nullTable)
    x <- segments[matchable, .(matchingSegment = unlist(matchingSegmentsTmp)), by = .(segmentName)]
    if(nrow(x) == 0) return(nullTable)    
    x[, i := 1:.N]
    x[, .(
        segment1 = min(segmentName, matchingSegment), # this paradigm for row-by-row reordering of two character columns is much faster
        segment2 = max(segmentName, matchingSegment) 
    ), keyby = .(i)][,
        .N,
        keyby = .(segment1, segment2) # two segments found to share at least one junction in common, thus are overlapping
    ]
}
analyzeSegmentsNetwork <- function(segments, segmentMatches){ # use igraph to identify and examine clusters of overlapping segments
    g <- graph_from_data_frame(segmentMatches, directed = FALSE)
    cmp <- components(g)
    matchable <- getMatchableSegments(segments)
    clusters <- merge( # scan the network of segment matches for clusters of molecules sharing junctions
        data.table(segmentName = names(cmp$membership), cluster = cmp$membership),
        segments[matchable, segmentLength, by = .(segmentName)],
        by = "segmentName",
        all.x = TRUE
    )
    setkey(clusters, cluster, segmentLength)
    clusters[, { # assign a reference segment for each cluster as the longest segment in the group
        .(
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
getSegmentPairs <- function(segmentGroups){ # work within the established segment groups to get potential overlap/identity pairs
    x <- segmentGroups[, data.table(t(combn(segmentName, 2))), by = .(refSegment)] 
    setnames(x, c("refSegment","segment1","segment2"))
    x
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

# examine segment pairs for duplex and identity relationships
getOuterEndpointMatch <- function(outerNodes1, outerNodes2){ # these junction nodes are original nodes, use fuzzy matching
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
scoreSegmentMatches <- function(segments, segmentPairs){
    do.call(rbind, mclapply(1:nrow(segmentPairs), function(i){
        segmentName1 <- segmentPairs[i, segment1]
        segmentName2 <- segmentPairs[i, segment2]
        segment1 <- segments[segmentName1]
        segment2 <- segments[segmentName2]
        data.table(
            refSegmentName = segmentPairs[i, refSegment],
            segmentName1 = segmentName1,
            segmentName2 = segmentName2,
            nodePathMatchType = getNodePathMatch(segment1[, unlist(nodePath)], segment2[, unlist(nodePath)]),
            outerEndpointMatchType = getOuterEndpointMatch(segment1[, c(outerNode1, outerNode2)], segment2[, c(outerNode1, outerNode2)]) 
        )
    }, mc.cores = env$N_CPU))  
}
