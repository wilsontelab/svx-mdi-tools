# -------------------------------------------------------------------------------------
# purge duplicate reads from edges
# -------------------------------------------------------------------------------------
# goal is to identify duplex reads to avoid having false elevation of the count of independent source molecules
# this is a read (not segment)-level assessment
# once done, every remaining read is taken as a true, independent source molecule (even if a chimeric one)
# (note that for amplification-free nanopore sequencing, reads from the same strand or different channels must come from independent molecules)
# -------------------------------------------------------------------------------------
# use the nanopore channel and the two outermost nodes to match reads to each other
# there is a rare but finite possibility that this process will purge independent molecules with the same random endpoints
# -------------------------------------------------------------------------------------
# after identifying duplicated duplex reads
#   entirely remove the duplicate read from edges (all junctions, all alignments)
#   update nStrands to 2 on the kept remaining reads
# -------------------------------------------------------------------------------------
# KNOWN LIMITATION - nanopore reads can often be partial spans, for reasons that
# might include strand breaks, low quality stretches, etc. This appears to happen frequently
# on paired duplex reads, e.g.:
#       ----------------|------------>
#       ~~~~<-----------|-------------
# where ~~~~ denotes bases presumed to be in the parent duplex that were not sequenced.
# Some insight into this situation can be inferred from having one shared endpoint, but even
# this is not wholly reliable. As a consequence, some true duplex reads of chimeric 
# molecules will persist and fail filtering, which demands that all fusable junctions 
# must have at least two unequivocally independent detections, either on the same strand
# or from different channels/nanopores.
# -------------------------------------------------------------------------------------

# concatenate nodes and edges into a single row per read (not segment)
collapseReads <- function(edges){
    message("collapsing edges into contiguous reads for duplex matching")
    reads <- edges[, {
        channel    <- channel[1]
        outerNode1 <- node1[1]
        outerNode2 <- node2[.N]
        isCanonical <- isCanonicalStrand(c(outerNode1, outerNode2))
        cOuterNode1 <- if(isCanonical) outerNode1 else -outerNode2
        cOuterNode2 <- if(isCanonical) outerNode2 else -outerNode1
        .(
            channel     = channel,
            outerNode1  = outerNode1, # original nodes of outer molecule endpoints, on strand(s) as sequenced
            outerNode2  = outerNode2,
            isCanonical = isCanonical, # strand orientation for reads is determined by the outermost alignments
            cOuterNode1 = cOuterNode1,
            cOuterNode2 = cOuterNode2,
            icPathKey   = paste(channel, cOuterNode1, cOuterNode2, sep = ":") 
        )        
    }, by = .(readI)] # a read-level assessment
    setkey(reads, icPathKey)
    reads
}

# examine duplex read relationships, i.e., same outer endpoints on opposite strands 
# identifies the minimal set of reads consistent with independent source molecules (not sequencing of the same molecule on both strands)
findDuplexReads <- function(reads){
    message("identifying apparent duplex reads")
    reads[, {
        nCanonical    <- sum( isCanonical)
        nNonCanonical <- sum(!isCanonical)
        if(nCanonical == 0 || nNonCanonical == 0) .(
            readI = readI,
            isDuplex = FALSE,
            retained = TRUE
        ) else .(
            readI = readI,
            isDuplex = TRUE, # keep all reads on the most abundant strand when duplex detections are found
            retained = isCanonical == (nCanonical >= nNonCanonical)
        )                 
    }, by = .(icPathKey)][, .SD, .SDcols = c("readI","isDuplex","retained")]
}
adjustEdgesForDuplex <- function(edges, duplexStatus){
    message("purging duplex reads")
    edges <- merge(
        edges,
        duplexStatus,
        by = "readI",
        all.x = TRUE
    )
    edges <- edges[retained == TRUE] # drops or keeps entire reads
    edges[isDuplex == TRUE, nStrands := 2]
    edges[, c("isDuplex","retained") := NULL]
    setkey(edges, readI, blockN, edgeN)
    edges
}
