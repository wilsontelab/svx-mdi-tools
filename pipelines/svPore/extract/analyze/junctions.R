# junction edges are pairs of nodes flanking a single junction
# use to assign reference junction designations as the most frequently used node pair among an adjacency cluster 
# as as first quick check for potential segment matches, i.e., if !any(junctionMatch), segmentMatch is not needed for a pair of segments
getJunctionsToMatch <- function(edges){
    if(sum(edges$keptJunction) == 0) stop("Zero (0) junctions were kept; sample has no apparent high-quality SV junctions")
    x <- edges[keptJunction == TRUE]
    x[, chromPairKey := paste(cChromIndex1, cChromIndex2, sep = ":")] # for more efficient adjacency searching
    x
}
getJunctionMatches <- function(junctions, junction_, junctionKey_){ # check a query edge against all other edges
    I <- junctions[, which(junctionKey != junctionKey_)]
    if(length(I) == 0) return(data.table(query = character(), match = character()))
    junctions[I,
        .( 
            match = junction_ %is_fuzzy_match% c(cNode1, cNode2),
            junctionKey = junctionKey
        ), 
        by = .(i)
    ][
        match == TRUE, 
        .(
            query = junctionKey_,
            match = junctionKey
        )
    ]
}
findMatchingJunctions <- function(junctionsToMatch){ # top-level function to perform all required junction match queries with parallelization
    x <- junctionsToMatch[,
        {
            junctions <- .SD[, .(cNode1, cNode2, junctionKey)]
            junctions[, i := 1:.N]
            do.call(rbind, mclapply(1:.N, function(i){
                .SD[i, getJunctionMatches(junctions, c(cNode1, cNode2), junctionKey)]
            }, mc.cores = env$N_CPU))
        }, 
        by = .(chromPairKey) # restrict search space to chrom pairs for speed
    ]
    if(nrow(x) == 0) return(data.table(junction1 = character(), junction2 = character()))
    x[, {
        x <- sort(c(query[1], match[1])) # order the matching junction edges (not the nodes within them)
        .(
            junction1 = x[1],
            junction2 = x[2] 
        )
    }, by = .(query, match)]
}
analyzeJunctionNetwork <- function(junctionMatches, junctionCounts){ # use igraph to identify and examine clusters of adjacent nodes
    if(nrow(junctionMatches) == 0) return(data.table( # handle case with no adjacencies
        junctionKey = character(), 
        isReferenceJunction = logical(), 
        rcNode1 = integer(), 
        rcNode2 = integer(),
        refJunctionKey = character()
    ))
    g <- graph_from_data_frame(junctionMatches[, .N, by = .(junction1, junction2)], directed = FALSE)
    cmp <- components(g)
    clusters <- merge( # scan the network of adjacencies for clusters of related junctions
        data.table(junctionKey = names(cmp$membership), cluster = cmp$membership),
        junctionCounts,
        by = "junctionKey",
        all.x = TRUE
    )
    setkey(clusters, cluster, nInstances)
    clusters[, { # assign a reference junction for each cluster as the the most frequent junction edge
        refJunctionNodes <- as.integer(strsplit(junctionKey[.N], ":")[[1]])
        .(
            junctionKey = junctionKey,
            isReferenceJunction = junctionKey == junctionKey[.N],        
            rcNode1 = refJunctionNodes[1], # "r" prefix means "reference", "c" again means "canonical" orientation
            rcNode2 = refJunctionNodes[2],
            refJunctionKey = paste(refJunctionNodes[1], refJunctionNodes[2], sep = ":")
        )
    }, by = .(cluster)][,
        .SD,
        .SDcols = c("junctionKey","isReferenceJunction","rcNode1","rcNode2","refJunctionKey")
    ]
}
finalizeNodeClustering <- function(edges, junctionClusters){ # assemble the final information for segment matching
    edges <- merge(edges, junctionClusters, by = "junctionKey", all.x = TRUE)
    edges[is.na(refJunctionKey), ":="( # fill in single-member junction clusters, i.e. with no adjancies (although junctionKey may have multiple segments)
        isReferenceJunction = TRUE, 
        rcNode1 = cNode1, 
        rcNode2 = cNode2,
        refJunctionKey = paste(cNode1, cNode2, sep = ":")
    )]
    edges[, ":="(
        rNode1 = ifelse(isCanonical, rcNode1, -rcNode2), # undo the canonical node flip
        rNode2 = ifelse(isCanonical, rcNode2, -rcNode1)
    )]
    edges[keptJunction == TRUE, ":="( # prepare for segment group by junction
        matchingSegments = list(unique(segmentName)),
        nInstances = .N
    ), by = .(refJunctionKey)]
    edges
}
