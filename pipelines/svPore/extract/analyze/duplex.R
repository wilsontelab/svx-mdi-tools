# -------------------------------------------------------------------------------------
# purge duplicate reads from edges
# -------------------------------------------------------------------------------------
# the goal is to identify duplex reads to avoid having false elevation of the count of independent source molecules
# this is a read (not segment)-level assessment
# once done, every remaining read is taken as a true source molecule (even if a chimeric one)
# (note that for amplification-free nanopore sequencing, reads on the same strand must come from independent molecules)
# -------------------------------------------------------------------------------------
# when defining the path used for read matching
#   use junction rNodes for more sensitive matching (fallback to nodes for adapter and low-quality junctions)
#   do NOT use junctions that failed bandwidth (they likely aren't real and may vary between duplex reads)
#   DO use junctions with adapters (has little impact, adapter chimeras will ~never be duplex)
#   DO use low quality junctions
# -------------------------------------------------------------------------------------
# after identifying duplicate reads
#   entirely remove them from edges (all junctions, all alignments)
# update nStrands to 2 on remaining reads when dropping duplex partners to a read
#   $ nStrands             : int  1 1 1 1 1 1 1 1 1 1 ...
# count the number of instances, i.e., unique source molecules, of each refJunctionKey 
# -------------------------------------------------------------------------------------

# concatenate nodes and edges into a single row per read (not segment)
parseReadNodePath <- function(real, node1, node2, rNode1, rNode2){
    rNode1 <- ifelse(is.na(rNode1), node1, rNode1) # fallback for real but unmatchable junctions that weren't scored by findMatchingJunctions
    rNode2 <- ifelse(is.na(rNode2), node2, rNode2)
    list(c(rbind(rNode1[real], rNode2[real])))
}
collapseReads <- function(edges){
    message("collapsing edges into contiguous reads for duplex matching")
    isReal <- getRealJunctions(edges)
    reads <- edges[, {
        .(
            outerNode1 = node1[1], # original nodes of outer molecule endpoints, on strand(s) as sequenced
            rNodePath = parseReadNodePath(isReal[.I], node1, node2, rNode1, rNode2), # reference junction node sequence, on strand(s) as sequenced
            outerNode2 = node2[.N]
        )
    }, by = .(qName)] # a read-level assessment
    reads[, ":="(i = 1:.N)]
    reads[, ":="( # strand orientation for reads is determined by the outermost alignments
        isCanonical = isCanonicalStrand(c(outerNode1, outerNode2))
    ), by = .(i)]
    reads[, ":="(
        cOuterNode1 = if(isCanonical) outerNode1 else -outerNode2,
        rcNodePath  = if(isCanonical) rNodePath  else list(-rev(unlist(rNodePath))),
        cOuterNode2 = if(isCanonical) outerNode2 else -outerNode1
    ), by = .(i)]
    reads[, ":="( # use exact matching on outer nodes for duplex purging
        pathKey = paste(cOuterNode1, unlist(rcNodePath), cOuterNode2, sep = ":", collapse = ":")
    ), by = .(i)]
    setkey(reads, pathKey)
    reads
}

# examine duplex read relationships, i.e., same node path on opposite strands, including outer endpoints
# identifies the minimal set of reads consistent with independent source molecules (not sequencing of the same molecule on both strands)
findDuplexReads <- function(reads){
    message("identifying apparent duplex reads")
    reads[, {
        nCanonical    <- sum( isCanonical)
        nNonCanonical <- sum(!isCanonical)
        if(nCanonical == 0 || nNonCanonical == 0) .(
            qName = qName,
            isDuplex = FALSE,
            retained = TRUE
        ) else .(
            qName = qName,
            isDuplex = TRUE, # keep all reads on the most abundant strand when duplex detections are found
            retained = isCanonical == (nCanonical >= nNonCanonical)
        )                 
    }, by = .(pathKey)]
}

adjustEdgesForDuplex <- function(edges, duplexStatus){
    message("purging duplex reads")
    edges <- merge(
        edges,
        duplexStatus,
        by = "qName",
        all.x = TRUE
    )
    edges <- edges[retained == TRUE] # drops or keeps entire reads
    edges[isDuplex == TRUE, nStrands := 2] 
    isJunction <- getJunctionEdges(edges)
    edges[isJunction, nJunctionInstances := .N, by = .(refJunctionKey)] # this count now represents unique source molecules
    edges
}
