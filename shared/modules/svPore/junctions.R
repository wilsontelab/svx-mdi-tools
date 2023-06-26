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
getJunctionsToMatch <- function(edges){
    if(edges[, sum(matchable, na.rm = TRUE) == 0]) stop("Edge list has no apparent high-quality SV junctions")
    x <- edges[matchable == TRUE]
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
        isIndexJunctionKey = logical(),
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
            junctionKey     = junctionKey,  # thus, output has one row per matched junctionKey, >=2 rows per fuzzy cluster
            isIndexJunctionKey = c(rep(FALSE, .N - 1), TRUE), # multiple junctions below may match this junctionKey
            nCanonical      = sum(nCanonical), # updated (potentially increased) counts AFTER considering window adjacencies
            nNonCanonical   = sum(nNonCanonical),
            icRefPos1       = indexData[1], # "i" prefix means "index junction of clusterN", "c" again means "canonical" orientation
            icRefPos2       = indexData[2],
            iInsertSize     = indexData[3]
        )
    }, by = .(clusterN)][,
        .SD,
        .SDcols = c("junctionKey","clusterN","isIndexJunctionKey","nCanonical","nNonCanonical","icRefPos1","icRefPos2","iInsertSize")
    ]
}
finalizeJunctionClustering <- function(edges, junctionHardCounts, junctionClusters){ # assemble the final information for read matching
    message("assigning junctions to clusters")

    # update edges within fuzzy-matched junction clusters
    #   as intended, this can merge a cluster identity into unmatchable junctions with the same junctionKey
    #   thereby updgrading them to "clustered" albeit "unmatchable"
    #   by virtue of sharing a junctionKey with a high quality junction
    edges <- merge(edges, junctionClusters, by = "junctionKey", all.x = TRUE)

    # initialize single-junctionKey clusters, i.e., with no fuzzy matches
    #   although the single junctionKey may hard-match multiple junction instances
    maxClusterN <- edges[, max(clusterN, na.rm = TRUE)]
    I <- edges[, matchable & is.na(clusterN)]
    clusterNs <- edges[I, .N, by = .(junctionKey)]
    clusterNs[, clusterN := 1:.N + maxClusterN]
    setkey(clusterNs, junctionKey) 
    edges[I, ":="( 
        clusterN = clusterNs[junctionKey, clusterN], # would it be faster to add the keyed items outside of groupby?  
        isIndexJunctionKey = TRUE,         
        nCanonical    = junctionHardCounts[junctionKey, nCanonical], 
        nNonCanonical = junctionHardCounts[junctionKey, nNonCanonical],  
        icRefPos1   = cRefPos1,
        icRefPos2   = cRefPos2,
        iInsertSize = insertSize        
    )]

    # set a convenience flag for tracking junction clusters along with matchable
    # as above, some "clustered" juntions may not have been "matchable"
    isJunction <- getJunctionEdges(edges)
    edges[isJunction, clustered := !is.na(clusterN)]

    # undo the canonical node flip
    edges[clustered == TRUE, ":="(
        iRefPos1 = ifelse(isCanonical, icRefPos1, icRefPos2), 
        iRefPos2 = ifelse(isCanonical, icRefPos2, icRefPos1)
    )]
    setkey(edges, sample, readI, blockN, edgeN)
    edges
}
collapseJunctionClusters <- function(edges){
    message("aggregating junctions by cluster group")

    # aggregate junctions to clusters
    junctions <- edges[clustered == TRUE, {
        i <- which(isIndexJunctionKey)[1] # select an index junction from those matching the indexJunctionKey
        samples <- unique(sample)
        .(
            edgeType        = edgeType[i], # these values are all the same for all instances of indexJunctionKey 
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
            mapQ            = max(mapQ), # these values aggregate over all junction in the cluster, even the fuzzy matched ones
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

    # stratify nInstances by sample
    junctionsBySample <- edges[clustered == TRUE, .(nInstances = .N), by = .(clusterN, sample)]
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
