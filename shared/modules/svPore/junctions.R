#-------------------------------------------------------------------------------------
# junction edges are pairs of nodes flanking a single junction
#-------------------------------------------------------------------------------------

# drop reads with no matchable junctions, they can never support an SV call
dropReadsWithNoJunctions <- function(edges, matchFilterFn){
    message("dropping reads with no usable junctions")
    matchable <- matchFilterFn(edges)
    I <- edges[, rep(any(matchable[.I]), .N), by = .(readI)][[2]] # keep or drop all edges from the read
    edges[I]
}

# perform fuzzy single-junction matching between reads
# used to assign reference junction designations as the most frequently used node pair among an adjacency cluster 
getJunctionsToMatch <- function(edges, matchFilterFn){
    matchable <- matchFilterFn(edges)
    if(nrow(edges) == 0 || sum(matchable) == 0) stop("Edge list has no apparent high-quality SV junctions")
    x <- edges[matchable]
    x[, searchKey := paste(cChromIndex1, cChromIndex2, cStrand1, cStrand2, sep = ":")] # for more efficient adjacency searching
    setkey(x, searchKey)
    x
}
getJunctionHardCounts <- function(junctionsToMatch){
    x <- junctionsToMatch[, .(
        nInstances      = .N,
        nCanonical      = sum( isCanonical), # BEFORE considering window adjancies, i.e., based on hard (exact) junction matching only
        nNonCanonical   = sum(!isCanonical)  # used to pick the index junction below
    ), by = .(junctionKey)]    
    setkey(x, junctionKey)
    x
}
findMatchingJunctions <- function(junctionsToMatch){ # top-level function to perform all required junction match queries with parallelization
    message("matching junctions between reads, with indel bandwidth tolerance")
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
        clusterN = integer(),        
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
        data.table(junctionKey = names(cmp$membership), clusterN = cmp$membership),
        junctionHardCounts,
        by = "junctionKey", 
        all.x = TRUE
    )
    setkey(clusters, clusterN, nInstances)
    clusters[, { # assign an index junction for each cluster as the the most frequently encountered junction edge
        indexData <- as.integer( strsplit(junctionKey[.N], ":")[[1]][5:7] )
        .(
            junctionKey     = junctionKey, 
            isIndexJunction = c(rep(FALSE, .N - 1), TRUE),
            nCanonical      = sum(nCanonical), # AFTER considering window adjacencies
            nNonCanonical   = sum(nNonCanonical),
            icRefPos1       = indexData[1], # "i" prefix means "index junction of clusterN", "c" again means "canonical" orientation
            icRefPos2       = indexData[2],
            iInsertSize     = indexData[3]
        )
    }, by = .(clusterN)][,
        .SD,
        .SDcols = c("junctionKey","clusterN","isIndexJunction","nCanonical","nNonCanonical","icRefPos1","icRefPos2","iInsertSize")
    ]
}
finalizeJunctionClustering <- function(edges, junctionHardCounts, junctionClusters){ # assemble the final information for read matching
    message("assigning junctions to clusters")

    # update edges within junction adjacency clusters
    # this may inappropriately merge a cluster into unmatchable junctions with the same key (see next step)
    edges <- merge(edges, junctionClusters, by = "junctionKey", all.x = TRUE)
    isJunction <- getJunctionEdges(edges)
    matchable <- getMatchableJunctions(edges)

    # override clusters for unmatchable junctions
    edges[isJunction & !matchable, ":="( 
        clusterN        = NA_integer_,   
        isIndexJunction = NA,
        nCanonical      = NA_integer_, 
        nNonCanonical   = NA_integer_,                 
        icRefPos1       = NA_integer_,
        icRefPos2       = NA_integer_,
        iInsertSize     = NA_integer_
    )]

    # initialize single-member junction clusters, i.e., with no adjacencies (although the single junctionKey may hard-match multiple SV instances)
    clusterN_ <- edges[, max(clusterN, na.rm = TRUE)]
    getClusterN <- function(){
        clusterN_ <<- clusterN_ + 1
        clusterN_
    }
    edges[matchable & is.na(clusterN), ":="( 
        clusterN        = getClusterN(),     
        isIndexJunction = c(rep(FALSE, .N - 1), TRUE),
        nCanonical      = junctionHardCounts[junctionKey, nCanonical], 
        nNonCanonical   = junctionHardCounts[junctionKey, nNonCanonical],                 
        icRefPos1       = cRefPos1[.N],
        icRefPos2       = cRefPos2[.N],
        iInsertSize     = insertSize[.N]
    ), by = .(junctionKey)]
    edges[matchable, ":="(
        iRefPos1 = ifelse(isCanonical, icRefPos1, icRefPos2), # undo the canonical node flip
        iRefPos2 = ifelse(isCanonical, icRefPos2, icRefPos1)
    )]
    setkey(edges, sample, readI, blockN, edgeN)
    edges
}
collapseJunctionClusters <- function(edges){
    message("aggregating junctions by cluster group")
    matchable <- getMatchableJunctions(edges)
    junctions <- edges[matchable, {
        i <- which(isIndexJunction)
        samples <- unique(sample)
        .(
            edgeType        = edgeType[i],
            eventSize       = eventSize[i],
            node1           = node1[i],
            node2           = node2[i],
            cChromIndex1    = cChromIndex1[i],
            cChromIndex2    = cChromIndex2[i],
            cStrand1        = cStrand1[i],
            cStrand2        = cStrand2[i],
            cRefPos1        = cRefPos1[i],
            cRefPos2        = cRefPos2[i],
            insertSize      = insertSize[i],
            mapQ            = max(mapQ),
            gapCompressedIdentity = max(gapCompressedIdentity),
            baseQual        = max(baseQual),
            alnBaseQual     = max(alnBaseQual),
            alnSize         = max(alnSize),
            samples         = paste(sort(samples), collapse = ","),
            nSamples        = length(samples),
            nInstances      = .N, # usually, but not necessarily == nReads (fails if same junction found >1 time in a read)
            nCanonical      = sum(isCanonical),
            nNonCanonical   = sum(!isCanonical)
        )
    }, by = .(clusterN)]
    junctionsBySample <- edges[matchable, .(nInstances = .N), by = .(clusterN, sample)]
    junctionsBySample <- dcast(junctionsBySample,
        clusterN ~ sample,
        fun.aggregate = sum,
        value.var = "nInstances"
    )
    merge(
        junctions,
        junctionsBySample,
        by = "clusterN",
        all.x = TRUE
    )
}
