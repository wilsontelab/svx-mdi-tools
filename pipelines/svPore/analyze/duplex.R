# -------------------------------------------------------------------------------------
# purge duplicate reads from edges
# -------------------------------------------------------------------------------------
# goal is to identify duplex reads to avoid having false elevation of the count of independent source molecules
# this is a read (not segment)-level assessment
# once done, every remaining read is taken as a true, independent source molecule (even if a chimeric one)
# (note that for amplification-free nanopore sequencing, reads on the same strand must come from independent molecules)
# -------------------------------------------------------------------------------------
# when defining the path used for read matching
#   use index junctions for more sensitive matching (fallback to nodes for adapter and low-quality junctions)
#   do NOT use junctions that failed bandwidth (they likely aren't real and may vary between duplex reads)
#   DO use junctions with adapters (has little impact, adapter chimeras will ~never be duplex)
#   DO use low quality junctions
# -------------------------------------------------------------------------------------
# after identifying duplicated duplex reads
#   entirely remove them from edges (all junctions, all alignments)
# update nStrands to 2 on remaining reads when dropping duplex partners to a read
# count the number of instances, i.e., remaining unique source molecules, of each indexJunctionKey 
# -------------------------------------------------------------------------------------
# KNOWN LIMITATION - unfortunately, nanopore reads can often be partial spans, for reasons that
# might include strand breaks, low quality stretches, etc. This appears to happen frequently
# on paired duplex reads, e.g.:
#       ----------------|------------>
#       ~~~~<-----------|-------------
# where ~~~~ denotes bases presumed to be in the parent duplex that were not sequenced.
# Some insight into this situation can be inferred from having one shared endpoint, but even
# this is probably not wholly reliable. As a consequence, some true duplex reads of chimeric 
# molecules will persist and fail filtering here, which demands that all fusable junctions 
# must have at least two detections _on the same strand_. Alternatively, could track and
# require that two different physical pores, i.e., nanopore channels, reported the same junction;
# that information is available in sequencing_summary.txt.
# -------------------------------------------------------------------------------------

# concatenate nodes and edges into a single row per read (not segment)
collapseReads <- function(edges, chromSizes){
    message("collapsing edges into contiguous reads for duplex matching")
    isReal <- getRealJunctions(edges)
    reads <- edges[, {
        .(
            outerNode1 = node1[1], # original nodes of outer molecule endpoints, on strand(s) as sequenced
            iReadPath = list(getIndexedReadPath(chromSizes, .SD, isReal[.I])), # index readPath, on strand(s) as sequenced
            outerNode2 = node2[.N],
            isCanonical = isCanonicalStrand(c(node1[1], node2[.N])) # strand orientation for reads is determined by the outermost alignments
        )
    }, by = .(qName)] # a read-level assessment
    reads[, ":="(
        cOuterNode1 = if(isCanonical) outerNode1 else -outerNode2,
        icReadPath  = list(getCanonicalReadPath(iReadPath[[1]], isCanonical)),
        cOuterNode2 = if(isCanonical) outerNode2 else -outerNode1
    ), by = .(qName)]
    reads[, ":="(
        icPathKey = getCanonicalPathKey(cOuterNode1, icReadPath[[1]], cOuterNode2) 
    ), by = .(qName)]
    setkey(reads, icPathKey)
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
    }, by = .(icPathKey)]
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
    edges[isJunction, ":="(
        nJunctionInstances = .N,        
        nCanonical    = sum( isCanonical), # recalculate canonicity for use by getFusableJunctions
        nNonCanonical = sum(!isCanonical)
    ), by = .(indexJunctionKey)] # this count now ~represents unique source molecules, but note LIMITATION above!
    setkey(edges, qName, blockN, edgeN)
    edges
}
