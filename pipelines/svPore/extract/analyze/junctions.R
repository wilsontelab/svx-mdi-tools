#-------------------------------------------------------------------------------------
# junction edges are pairs of nodes flanking a single junction
#-------------------------------------------------------------------------------------

# drop reads with no matchable junctions, they can never support an SV call
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
    x[, searchKey := paste(cChromIndex1, cChromIndex2, cStrand1, cStrand2, sep = ":")] # for more efficient adjacency searching
    setkey(x, searchKey)
    x
}
getJunctionHardCounts <- function(junctionsToMatch){
    x <- junctionsToMatch[, .(
        nInstances      = .N,
        nCanonical      = sum( isCanonical), # BEFORE considering window adjancies, i.e., based on hard (exact) junction matching only
        nNonCanonical   = sum(!isCanonical)  # used to pick the reference junction below
    ), by = .(junctionKey)]    
    setkey(x, junctionKey)
    x
}
findMatchingJunctions <- function(junctionsToMatch){ # top-level function to perform all required junction match queries with parallelization
    message("matching junctions between reads, with junction indel tolerance")
    x <- junctionsToMatch[,
        {
            searchJunctions <- .SD[, .(junctionKey, cStrand1, cStrand2, cRefPos1, cRefPos2, insertSize)]
            do.call(rbind, mclapply(1:.N, function(j){
                junctionKey_ <- searchJunctions[j, junctionKey]
                acRefPos1_   <- searchJunctions[j, cRefPos1 + cStrand1 * insertSize / 2] # "a" means "adjusted for insertSize"
                acRefPos2_   <- searchJunctions[j, cRefPos2 - cStrand2 * insertSize / 2]
                searchJunctions[, {
                    acRefPos1 <- cRefPos1 + cStrand1 * insertSize / 2
                    acRefPos2 <- cRefPos2 - cStrand2 * insertSize / 2
                    .(
                        queryJxnKey = junctionKey_,
                        targetJxnKey = junctionKey,
                        distance = sqrt((acRefPos1_ - acRefPos1)**2 + (acRefPos2_ - acRefPos2)**2)
                    )
                }]
            }, mc.cores = env$N_CPU))
        }, 
        by = .(searchKey) # restrict search space to chrom pairs for speed
    ]
    if(nrow(x) == 0) return(data.table(junctionKey1 = character(), junctionKey2 = character()))
    x[, k := 1:.N]
    x[queryJxnKey != targetJxnKey & distance <= env$JUNCTION_BANDWIDTH, .( # order the matching junction edges (not the nodes within them)
        junctionKey1 = min(queryJxnKey, targetJxnKey), # this paradigm for row-by-row reordering of two character columns is much faster
        junctionKey2 = max(queryJxnKey, targetJxnKey)
    ), keyby = .(k)][,
        .N,
        keyby = .(junctionKey1, junctionKey2) # two junctions that fuzzy-matched across at least two reads that need to be associated with each other
    ]
}
analyzeJunctionNetwork <- function(junctionMatches, junctionHardCounts){ # use igraph to identify and examine clusters of adjacent nodes
    if(nrow(junctionMatches) == 0) return(data.table( # handle case with no fuzz-match junction adjacencies
        junctionKey = character(), 
        indexJunctionKey = character(),        
        isIndexJunction = logical(), 
        nCanonical = integer(),
        nNonCanonical = integer(),        
        icRefPos1 = integer(), 
        icRefPos2 = integer(),
        iInsertSize = integer()
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
    clusters[, { # assign a reference junction for each cluster as the the most frequently encountered junction edge
        indexData <- as.integer( strsplit(junctionKey[.N], ":")[[1]][5:7] )
        .(
            junctionKey = junctionKey,
            indexJunctionKey = junctionKey[.N],             
            isIndexJunction = junctionKey == junctionKey[.N],   
            nCanonical = sum(nCanonical), # AFTER considering window adjancies
            nNonCanonical = sum(nNonCanonical),
            icRefPos1   = indexData[1], # "i" prefix means "index junction of cluster", "c" again means "canonical" orientation
            icRefPos2   = indexData[2],
            iInsertSize = indexData[3]
        )
    }, by = .(cluster)][,
        .SD,
        .SDcols = c("junctionKey","indexJunctionKey","isIndexJunction","nCanonical","nNonCanonical","icRefPos1","icRefPos2","iInsertSize")
    ]
}
finalizeJunctionClustering <- function(edges, junctionHardCounts, junctionClusters){ # assemble the final information for read matching
    message("assigning reference junctions")

    # update edges within adjacency clusters
    edges <- merge(edges, junctionClusters, by = "junctionKey", all.x = TRUE)
    matchable <- getMatchableJunctions(edges)

    # update single-member junction clusters, i.e. with no adjacencies (although junctionKey may hard-match multiple SV instances)
    edges[matchable & is.na(indexJunctionKey), ":="( 
        indexJunctionKey = junctionKey,             
        isIndexJunction = TRUE,   
        nCanonical = NA_integer_, 
        nNonCanonical = NA_integer_,                 
        icRefPos1   = cRefPos1,
        icRefPos2   = cRefPos2,
        iInsertSize = insertSize
    )]
    edges[matchable & is.na(nCanonical), ":="( 
        nCanonical    = junctionHardCounts[indexJunctionKey, nCanonical],
        nNonCanonical = junctionHardCounts[indexJunctionKey, nNonCanonical]
    )]
    edges[matchable, ":="(
        iRefPos1 = ifelse(isCanonical, icRefPos1, icRefPos2), # undo the canonical node flip
        iRefPos2 = ifelse(isCanonical, icRefPos2, icRefPos1)
    )]
    setkey(edges, qName, blockN, edgeN)
    edges
}

# extract and analyze all singleton junctions
pullSingletonJunctions <- function(edges){
    singletons <- getSingletonJunctions(edges)
    
}
