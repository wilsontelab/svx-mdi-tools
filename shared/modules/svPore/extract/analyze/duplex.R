# -------------------------------------------------------------------------------------
# purge duplicate duplex read pairs from edges
# -------------------------------------------------------------------------------------
# from explorations of `dorado duplex` output summarized in part here:
#   https://github.com/nanoporetech/dorado/issues/443
# the following summarizes the disposition of simplex and duplex reads through basecalling and alignment
#   - most basecalled reads are simplex/unsplit 
#   - split reads do not necessarily proceed to duplexing (i.e., splitting is based on adapters, not self-complementarity)
#   - foldback chimeras with adapters can be split and duplexed
#   - duplexes yielding two different unsplit reads can also be found and duplexed
#   - HOWEVER, unsplittable foldback chimeras lacking internal adapters are NOT found as duplexes and persist in Dorado output
# ------------------------------------------------------------------------------------
# THUS, extract/parse_node.pl:checkForDuplex() used minimap2 alignments to purge 
# those false inversions and flag the remaining reads as foldback==1 duplex==0
# to reflect that they were sequenced as duplexes but failed to undergo stereo basecalling
# foldback reads are single reads where the foldback redundancy and false inversion junction have been purged
# -------------------------------------------------------------------------------------
# but we aren't done yet!
# here we continue to identify duplex read pairs to avoid having false elevation of the count of independent source molecules
# this is a read (not segment)-level assessment
# once done, every non-redunant read is taken as a true, independent source molecule (even if a chimeric one)
# -------------------------------------------------------------------------------------
# logic is to use the nanopore channel and the two outermost nodes to fuzzy-match reads to each other
# there is a rare but finite possibility this process will flag independent molecules with the same random endpoints
# -------------------------------------------------------------------------------------
# after identifying duplex molecules as duplicated simplex reads
#   set column 'duplex2' to 0 for non-duplex reads, -1 for duplex reads, to match the usage in Dorado found the the 'duplex' column
#   all simplex reads identified by either dorado or us are retained, i.e., SVs are still duplicated in simplex reads
#   finally set duplexCluster and duplexCluster2 to a value that links simplex pairs to each other
# -------------------------------------------------------------------------------------
# KNOWN LIMITATION - nanopore reads can often be partial spans, for reasons that
# might include strand breaks, low quality stretches, etc. This can happen on paired duplex reads, e.g.:
#       ----------------|------------> read 1, channel 1
#       ~~~~<-----------|---------~~~~ read 2, channel 1
# where ~~~~ denotes bases presumed to be in the parent duplex that were not sequenced or basecalled
# As a consequence, some true duplex reads will persist and fail filtering even through here. 
# Later, collapseJunctionClusters only counts junctions as independent if they 
# derive from the same strand or different channels/nanopores.
# -------------------------------------------------------------------------------------

# lookup the read id(s) of the Dorado duplex reads derived from two simplex reads
assignDoradoDuplexClusters <- function(edges, duplexEdges){
    message("linking Dorado simplex read pairs in duplex reads")
    duplexClusters <- duplexEdges[, .(
        simplexQName = strsplit(qName, ";")[[1]]
    ), by = .(duplexCluster)]
    setkey(duplexClusters, simplexQName)
    edges[
        duplex == -1, 
        duplexCluster := paste(duplexClusters[qName, duplexCluster], collapse = ":"),
        by = .(qName)
    ]
    edges
}

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

    message("marking and linking duplex reads")
    setkey(clusters, clusterN, readI)    
    clusterKeys <- clusters[, .(clusterKey = paste(readI, collapse = ":")), by = .(clusterN)]
    clusterKeys[, clusterN_ := as.character(clusterN)]
    setkey(clusterKeys, clusterN_)
    clusters[, readI_ := as.character(readI)]
    setkey(clusters, readI_)
    edges[readI %in% clusters$readI, ":="(
        duplex2 = -1,
        duplexCluster2 = {
            readI_ <- as.character(readI)
            clusterN_ <- clusters[readI_, as.character(clusterN)]
            clusterKeys[clusterN_, clusterKey]
        }
    )]
    setkey(edges, readI, blockN, edgeN)
    edges
}
