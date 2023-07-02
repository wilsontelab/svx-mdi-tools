# -------------------------------------------------------------------------------------
# purge duplicate duplex read pairs from edges
# -------------------------------------------------------------------------------------
# goal is to identify duplex read pairs to avoid having false elevation of the count of independent source molecules
# this is a read (not segment)-level assessment
# once done, every remaining read is taken as a true, independent source molecule (even if a chimeric one)
# -------------------------------------------------------------------------------------
# this code is designed to find duplex read pairs that were called as separate reads
# upstream code in extract/parse_node.pl:checkForDuplex() filters against inline inversion foldbacks called as a single read
# -------------------------------------------------------------------------------------
# use the nanopore channel and the two outermost nodes to fuzzy-match reads to each other
# there is a rare but finite possibility that this process will purge independent molecules with the same random endpoints
# -------------------------------------------------------------------------------------
# after identifying duplicated duplex reads
#   entirely remove the duplicate read from edges (all junctions, all alignments)
#   update nStrands to 2 on the kept remaining reads
# -------------------------------------------------------------------------------------
# KNOWN LIMITATION - nanopore reads can often be partial spans, for reasons that
# might include strand breaks, low quality stretches, etc. This can happen on paired duplex reads, e.g.:
#       ----------------|------------>
#       ~~~~<-----------|-------------
# where ~~~~ denotes bases presumed to be in the parent duplex that were not sequenced.
# Some insight into this situation can be inferred from having one shared endpoint, but even
# this is not wholly reliable. As a consequence, some true duplex reads of chimeric 
# molecules will persist and fail filtering here. Later, collapseJunctionClusters only counts 
# junctions as independent if they derive from the same strand or different channels/nanopores.
# -------------------------------------------------------------------------------------

# concatenate edges into a single row per read (not segment)
collapseReads <- function(edges, chromSizes){
    message("collapsing edges into contiguous reads for duplex matching")
    reads <- edges[, {
        outerNode1  <- node1[1]
        outerNode2  <- node2[.N]
        isCanonical <- isCanonicalStrand(c(outerNode1, outerNode2))
        .(
            channel     = channel[1],
            cOuterNode1 = if(isCanonical) outerNode1 else -outerNode2,
            cOuterNode2 = if(isCanonical) outerNode2 else -outerNode1,
            nJunctions  = sum(edgeType != edgeTypes$ALIGNMENT) 
        )        
    }, by = .(readI)] # a read-level assessment
    reads <- cbind(
        reads,
        reads[, parseSignedNodes(chromSizes, cOuterNode1, 1, canonical = TRUE)], 
        reads[, parseSignedNodes(chromSizes, cOuterNode2, 2, canonical = TRUE)]
    )
    reads[, searchKey := paste(channel, cChromIndex1, cChromIndex2, cStrand1, cStrand2, sep = ":")] # duplex matching is channel-restricted and assisted by search space
    setkey(reads, searchKey)
    reads
}

# perform all required read match queries with parallelization
findMatchingReads <- function(reads){ 
    message("matching reads based on outer endpoints, with alignment tolerance")
    nullXTable <- data.table(
        queryReadI  = integer(),
        targetReadI = integer(),
        distance    = double()
    )
    x <- reads[,
        {
            if(.N == 1) nullXTable # most often, there is at most one duplex match within the read search space (unlike junctions)
            else if(.N == 2) data.table(
                queryReadI  = readI[1],
                targetReadI = readI[2],
                distance    = sqrt((cRefPos1[1] - cRefPos1[2])**2 + (cRefPos2[1] - cRefPos2[2])**2)
            ) else { 
                searchReads <- .SD[, .(readI, cRefPos1, cRefPos2)]
                do.call(rbind, mclapply(1:.N, function(j){
                    readI_     <- searchReads[j, readI]
                    cRefPos1_  <- searchReads[j, cRefPos1]
                    cRefPos2_  <- searchReads[j, cRefPos2]
                    searchReads[, {
                        .(
                            queryReadI  = readI_,
                            targetReadI = readI,
                            distance = sqrt((cRefPos1_ - cRefPos1)**2 + (cRefPos2_ - cRefPos2)**2)
                        )
                    }]
                }, mc.cores = env$N_CPU))                
            }
        }, 
        by = .(searchKey) # restrict search space to chrom pairs for speed
    ]
    nullTable <- data.table(readI1 = integer(), readI2 = integer())
    if(nrow(x) == 0) return(nullTable)
    x[, ":="(
        k = 1:.N,
        passed = queryReadI != targetReadI & distance <= env$JUNCTION_BANDWIDTH
    )]
    if(x[, sum(passed) == 0]) return(nullTable)
    x[passed == TRUE, .( # order the non-self matching reads
        readI1 = min(queryReadI, targetReadI), # this paradigm for row-by-row reordering of two character columns is much faster
        readI2 = max(queryReadI, targetReadI)
    ), keyby = .(k)][,
        .N,
        keyby = .(readI1, readI2) # two different reads that fuzzy-matched their outer endpoints from the same nanopore channel, presumably duplex pairs
    ]
}
analyzeReadsNetwork <- function(readMatches, reads, edges){ # use igraph to identify and examine clusters of adjacent reads
    message("analyzing the duplex reads network")
    if(nrow(readMatches) == 0) return(edges)
    g <- graph_from_data_frame(readMatches, directed = FALSE)
    cmp <- components(g)
    clusters <- merge( # scan the network of adjacencies for clusters of related reads
        data.table(readI = as.integer(names(cmp$membership)), clusterN = cmp$membership),
        reads,
        by = "readI", 
        all.x = TRUE
    )
    setkey(clusters, clusterN, nJunctions)

    message("purging duplex reads")
    indexReadIs  <- clusters[, readI[1],    by = .(clusterN)] # for each cluster, keep the read with the fewest junctions
    duplexReadIs <- clusters[, readI[2:.N], by = .(clusterN)]
    edges <- edges[!(readI %in% duplexReadIs)] # drop all non-index reads
    edges[readI %in% indexReadIs, nStrands := 2]
    setkey(edges, readI, blockN, edgeN)
    edges
}
