#-------------------------------------------------------------------------------------
# svPore analyze support functions
#-------------------------------------------------------------------------------------
# load input data
loadEdges <- function(type) {
    message(paste("loading", type, "edges"))
    fread(
        if(type == "sv") env$EDGES_SV_FILE else env$EDGES_TMP_FILE,
        col.names = edgesCols,
        colClasses = edgesColClasses,
        sep = "\t",
        quote = ""
    )
}
loadReads <- function(){
    message("loading read sequences")
    fread(
        env$SEQUENCES_FILE,
        col.names = readsCols,
        colClasses = readsColClasses,
        sep = "\t",
        quote = ""
    )
}

# expand integer nodes out to chrom/strand/pos
parseSignedWindow <- function(window, side) {
    strand <- ifelse(window > 0, "+", "-")
    window <- abs(window)
    chromIndex <- bitwShiftR(window, 24)
    dt <- data.table(
        chromIndex  = chromIndex,
        chrom       = unlist(revChromIndex[chromIndex]),
        windowIndex = bitwAnd(window, 2**24 - 1),
        strand      = strand
    )
    names(dt) <- paste0(names(dt), side)
    dt
}

# rectify a junction edge, i.e., a pair of nodes, to a consistently ordered string representation suitable for sorting and grouping
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0
}
sortCanonical <- function(nodes){
    isCanonical <- isCanonicalStrand(c(nodes[1], nodes[length(nodes)]))
    if(!isCanonical) nodes <- -rev(nodes)
    nodes
}
getPathSignature <- function(nodes, reverse = FALSE){
    isCanonical <- isCanonicalStrand(c(nodes[1], nodes[length(nodes)]))
    if(!isCanonical) nodes <- -rev(nodes)
    if(reverse) nodes <- -rev(nodes)
    paste0(":", paste(nodes, collapse = ":"), ":")
}

# functions that support:
#   fuzzy SV path matching, i.e., with tolerance for matching nodes in closely adjacent windows
#   partial SV path matching, i.e., that reports on molecule overlap/nesting as well as duplex/identity
pathTolerance <- 1
pathOffsets <- -pathTolerance:pathTolerance
`%is_path%` <- function(path1, path2){ # core function that checks if two equal-length node sequences are identical
    all(apply(sapply(pathOffsets, "+", path2) - path1, 1, function(v) min(abs(v))) == 0)
}
`%is_rev_path%` <- function(path1, path2){ # check if two paths are exactly reciprocal with respect to strand
    path1 %is_path% -rev(path2)
}

# junction edges are pairs of nodes flanking the junction itself
# used as first quick check for potential segment matches, i.e., if !any(getJunctionMatch), getSegmentMatch is not needed
junctionMatchTyes <- list(
    NONE            = 0,
    SAME_STRAND     = 1,
    OPPOSITE_STRAND = 2
)
getJunctionMatch <- function(edge1, edge2){ 
    if(edge1 %is_rev_path% edge2)  junctionMatchTyes$OPPOSITE_STRAND 
    else if(edge1 %is_path% edge2) junctionMatchTyes$SAME_STRAND
    else junctionMatchTyes$NONE
}
getJunctionMatches <- function(segmentName1, blockN1, edgeN1, edgeType_, chromPair_, edge1){ # check a query edge against all other edges, including itself as a sanity check
    x <- junctionsToMatch[
        edgeType == edgeType_ & # single junctions must have a common edgeType (whereas molecule mights have different composite pathTypes)
        chromPair == chromPair_ & # restrict the search space to the chromosome pair for speed
        segmentName != segmentName1 # don't match molecules to themselves        
    ]
    if(nrow(x) == 0) return(NULL)
    x[,
        .( matchType = getJunctionMatch(edge1, c(node1, node2)) ), 
        by = .(segmentName, blockN, edgeN)
    ][
        matchType > junctionMatchTyes$NONE, 
        .(
            segmentName = segmentName1,
            blockN = blockN1,
            edgeN = edgeN1,
            matchingSegments = list(segmentName)
        )
    ]
}

# segments are arbitrary length vectors of high quality nodes, including the outer molecule endpoints
# there is always an even number of nodes as passed to these functions (2 per junction, 2 outer alignment nodes)
pathMatchTypes <- list(
    NONE        = 0,
    IDENTICAL   = 1, # IDENTICAL, DUPLEX, and SAME_NODES demand an equal number of nodes in the two paths
    DUPLEX      = 2, # IDENTICAL and DUPLEX require matching outer endpoints, i.e., molecules are considered to be the same source DNA molecule
    SAME_NODES  = 3, # SAME_NODES, NESTED, and OVERLAP do NOT consider molecule outer endpoints, i.e., molecules are considered to be different source DNA molecules
    NESTED      = 4, # NESTED has one segment entirely contained in another
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
        if(x1[i] %is_path% x2[i] || x1r[i] %is_path% x2[i]) return(pathMatchTypes$OVERLAP)
    }
    pathMatchTypes$NONE
}
getSegmentMatch <- function(segment1, segment2){ # segment must be c(outerNode, junction, ..., outerNode)
    path1 <- segment1[2:(length(segment1) - 1)]
    path2 <- segment2[2:(length(segment2) - 1)]
    if(length(segment1) == length(segment2)){
        if(path1 %is_rev_path% path2 || path1 %is_path% path2) { # check without outer endpoints first, since most targets will not be a match
            if(segment1 %is_rev_path% segment2)  pathMatchTypes$DUPLEX # some duplication, but this is an infrequent code path
            else if(segment1 %is_path% segment2) pathMatchTypes$IDENTICAL            
            else pathMatchTypes$SAME_NODES
        } else if(length(path1) > 2 || length(path2) > 2) getOverlapMatch(path1, path2)
        else pathMatchTypes$NONE
    } else getOverlapMatch(path1, path2)
}
getSegmentMatches <- function(segmentName1, matchingSegments_, segment1){ # check a query edge against all other edges, including itself as a sanity check
    if(is.null(matchingSegments_)) return(NULL)
    segmentsToMatch[
        segmentName %in% matchingSegments_, # the list of segments identified above as having >=1 junction match to segment1
        .( matchType = getSegmentMatch(segment1, unlist(nodePath)) ), 
        by = .(segmentName)
    ][, .(
        segmentName = segmentName1,
        matchingSegments = list(paste(segmentName, matchType, sep = ":"))
    )]
}
