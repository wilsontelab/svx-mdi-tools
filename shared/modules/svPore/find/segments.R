#-------------------------------------------------------------------------------------
# segments are arbitrary length series of high quality nodes, including the outer molecule endpoints
# that describe a single, contiguous, confidently _non-chimeric_ portion of a read
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
splitReadsToSegments <- function(edges, fusableFn = NULL, msg = NULL){
    if(is.null(fusableFn)) fusableFn <- getFusableJunctions
    if(is.null(msg)) msg <- "splitting reads at low-quality, chimeric, or unconfirmed junctions"
    message(msg)    
    fusable <- fusableFn(edges)
    edges[, segmentN := getSplitSegmentN(.N, !fusable[.I]), by = .(sample, readI)]
    edges <- edges[!is.na(segmentN)] # drop the false/untrusted junctions
    setkey(edges, sample, readI, segmentN)

    # update blockN and edgeN within each split segment
    edges[, ":="(
        segmentName = paste(qName, segmentN, sep = "-"), # segmentName is unique even between samples since qName is a uuid
        blockN = blockN - min(blockN) + 1,
        edgeN  = edgeN  - min(edgeN)  + 1
    ), by = .(sample, readI, segmentN)]

    # remove split segments that do not contain any potential SV junctions
    # e.g., in A[T]A[T]A or A[T]AdA[T]A (where "d" was rejected but not split)
    # TODO: save these segments as alignments for haplotype assembly?
    isClustered <- "clusterN" %in% names(edges)
    clustered_ <- if(isClustered) edges$clustered else TRUE
    segmentsToKeep <- edges[clustered_, unique(segmentName)] # segments with at least one clustered junction
    edges <- edges[segmentName %in% segmentsToKeep]
    setkey(edges, segmentName, blockN, edgeN)

    # prepare for assembly segment matching
    if(isClustered) edges[clustered == TRUE, ":="( 
        matchingSegments = list(unique(segmentName))
    ), by = .(clusterN)]
    edges
}

# concatenate nodes and edges into a single row per segment
# every segment has at least one clustered junction with nInstances > 1
collapseSegments <- function(edges, chromSizes){
    message("collapsing edges into contiguous, trusted read segments")
    isAlignment <- getAlignmentEdges(edges)
    segments <- edges[, {
        alignment <- isAlignment[.I]
        clustered <- !is.na(clustered) & clustered
        .(
            sample = sample[1],
            readI = readI[1],
            segmentLength = qEnd[.N] - qStart[1],
            chroms = paste(unique(c(chrom1, chrom2)), collapse = ","),
            nJxns = sum(!alignment),
            nClusteredJxns = sum(clustered, na.rm = TRUE),
            pathType = paste0(ifelse(alignment | clustered, edgeType, tolower(edgeType)), collapse = ""),
            isClosedPath = chromIndex1[1] == chromIndex2[.N] && strand1[1] == strand2[.N],
            outerNode1 = node1[1], # original nodes of outer molecule endpoints, on strand(s) as sequenced
            iReadPath = list(getIndexedReadPath(chromSizes, .SD, clustered)), # index readPath, on strand(s) as sequenced
            outerNode2 = node2[.N],
            isCanonical = isCanonicalStrand(c(node1[1], node2[.N])), # strand orientation for segments is determined by the outermost alignments
            matchingSegmentsTmp = { # a matching segment has at least one fuzzy junction in common with the query segment
                x <- unique(unlist(matchingSegments))
                list(x[x != segmentName])
            }
        )
    }, by = .(segmentName)]
    segments[, ":="(
        nMatchingSegments = sapply(matchingSegmentsTmp, function(x) length(unlist(x))),
        cOuterNode1 = if(isCanonical) outerNode1 else -outerNode2,
        icReadPath  = list(getCanonicalReadPath(iReadPath[[1]], isCanonical)),
        cOuterNode2 = if(isCanonical) outerNode2 else -outerNode1
    ), by = .(segmentName)]
    segments[, ":="( # use exact matching on outer nodes for duplex purging
        icPathKey = getCanonicalPathKey(cOuterNode1, icReadPath[[1]], cOuterNode2) 
    ), by = .(segmentName)]
    segments
}

# assemble the networks of segments with shared junctions, i.e., overlapping segments on a haplotype
getNuclearSegments <- function(segments){
    message("removing segments confined to chrM (only nuclear genome SVs are reported for assembly)")
    x <- segments[chroms != "chrM"]
    setkey(x, segmentName)
    x
}
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
        data.table(    # these clusters do _not_ necessarily correspond to a single heplotype in a mosaic/diploid sample
            segmentName = names(cmp$membership), 
            clusterN = cmp$membership
        ),
        segments[matchable, .(
            nClusteredJxns,
            segmentLength, 
            pathType
        ), by = .(segmentName)],
        by = "segmentName",
        all.x = TRUE
    )
    setkey(clusters, clusterN, segmentLength, nClusteredJxns)
    clusters[, { # assign an index segment for each cluster as the longest segment in the group of segments with the most clustered junctions
        .(
            segmentName = segmentName,
            isIndexSegment = c(rep(FALSE, .N - 1), TRUE),
            pathType = pathType
        )
    }, by = .(clusterN)][,
        .SD,
        .SDcols = c("segmentName","clusterN","isIndexSegment","pathType") 
    ]
}
collapseSegmentClusters <- function(segmentEdges, segmentsNetwork){
    message("aggregating segments by cluster group")
    clusteredEdges <- segmentEdges[clustered == TRUE]
    setkey(clusteredEdges, segmentName)
    clusters <- segmentsNetwork[, {
        clusterSegmentNames <- unique(segmentName)
        junctionClusterNs <- clusteredEdges[clusterSegmentNames, unique(na.omit(clusterN))]   
        .(
            nSegments = .N,
            segmentNames = list(clusterSegmentNames),
            pathTypes = paste(sort(unique(pathType)), collapse = " "),
            nJunctionClusters = length(junctionClusterNs),
            junctionClusterNs = list(junctionClusterNs)     
        )        
    }, by = .(clusterN)]
    clusters
}
updateJunctionClustersForSegments <- function(junctionClusters, segmentClusters){
    message("updating junction clusters with segment cluster membership")    
    x <- segmentClusters[, .(junctionClusterN = unlist(junctionClusterNs)), by = .(clusterN)] # a junction cluster can be in only one segment cluster
    setnames(x, c("segmentClusterN", "clusterN"))
    merge(junctionClusters, x, by = "clusterN", all.x = TRUE)
}
