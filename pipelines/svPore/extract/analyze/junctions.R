#-------------------------------------------------------------------------------------
# junction edges are pairs of nodes flanking a single junction
#-------------------------------------------------------------------------------------

# drop reads with no matchable junctions
dropReadsWithNoJunctions <- function(edges){
    message("dropping reads with no matchable junctions")
    matchable <- getMatchableJunctions(edges)
    I <- edges[, rep(any(matchable[.I]), .N), by = .(qName)][[2]]
    edges[I]
}

# perform fuzzy single-junction matching between reads
# used to assign reference junction designations as the most frequently used node pair among an adjacency cluster 
getJunctionsToMatch <- function(edges){
    matchable <- getMatchableJunctions(edges)
    if(nrow(edges) == 0 || sum(matchable) == 0) stop("Sample has no apparent high-quality SV junctions")
    x <- edges[matchable]
    x[, chromPairKey := paste(cChromIndex1, cChromIndex2, sep = ":")] # for more efficient adjacency searching
    x
}
getJunctionHardCounts <- function(junctionsToMatch){
    x <- junctionsToMatch[, .(
        nInstances      = .N,
        nCanonical      = sum(isCanonical), # BEFORE considering window adjancies, i.e., based on hard junction matching only
        nNonCanonical   = sum(!isCanonical)
    ), by = .(junctionKey)]    
    setkey(x, junctionKey)
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
    message("matching junctions between reads, with adjacency tolerance")
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
    x[, i := 1:.N]
    x[, .( # order the matching junction edges (not the nodes within them)
        junction1 = min(query, match), # this paradigm for row-by-row reordering of two character columns is much faster
        junction2 = max(query, match) 
    ), keyby = .(i)][,
        .N,
        keyby = .(junction1, junction2) # two segments found to share at least one junction in common, thus are overlapping
    ]
}
analyzeJunctionNetwork <- function(junctionMatches, junctionHardCounts){ # use igraph to identify and examine clusters of adjacent nodes
    if(nrow(junctionMatches) == 0) return(data.table( # handle case with no adjacencies
        junctionKey = character(), 
        isReferenceJunction = logical(), 
        rcNode1 = integer(), 
        rcNode2 = integer(),
        refJunctionKey = character(),
        nCanonical = integer(),
        nNonCanonical = integer()
    ))
    g <- graph_from_data_frame(junctionMatches, directed = FALSE)
    cmp <- components(g)
    clusters <- merge( # scan the network of adjacencies for clusters of related junctions
        data.table(junctionKey = names(cmp$membership), cluster = cmp$membership),
        junctionHardCounts,
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
            refJunctionKey = paste(refJunctionNodes[1], refJunctionNodes[2], sep = ":"),
            nCanonical = sum(nCanonical), # AFTER considering window adjancies
            nNonCanonical = sum(nNonCanonical)
        )
    }, by = .(cluster)][,
        .SD,
        .SDcols = c("junctionKey","isReferenceJunction","rcNode1","rcNode2","refJunctionKey","nCanonical","nNonCanonical")
    ]
}
finalizeJunctionClustering <- function(edges, junctionHardCounts, junctionClusters){ # assemble the final information for read matching
    message("assigning reference junctions")
    edges <- merge(edges, junctionClusters, by = "junctionKey", all.x = TRUE)
    matchable <- getMatchableJunctions(edges)
    # fill in single-member junction clusters, i.e. with no adjancies (although junctionKey may hard-match multiple instances)
    edges[matchable & is.na(refJunctionKey), ":="( 
        isReferenceJunction = TRUE, 
        rcNode1 = cNode1, 
        rcNode2 = cNode2,
        refJunctionKey = paste(cNode1, cNode2, sep = ":"),
        nCanonical = NA_integer_,
        nNonCanonical = NA_integer_
    )]
    edges[matchable & is.na(nCanonical), ":="( 
        nCanonical    = junctionHardCounts[refJunctionKey, nCanonical],
        nNonCanonical = junctionHardCounts[refJunctionKey, nNonCanonical]
    )]
    edges[matchable, ":="(
        rNode1 = ifelse(isCanonical, rcNode1, -rcNode2), # undo the canonical node flip
        rNode2 = ifelse(isCanonical, rcNode2, -rcNode1)
    )]
    edges
}
