#-------------------------------------------------------------------------------------
# filters.R returns filters of edges and segments based on their assessed properties
# all functions return a logical vector of the same length as the provided data.table
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# edge filters
#-------------------------------------------------------------------------------------

# all alignment or junction edges
getAlignmentEdges <- function(edges) edges[, edgeType == edgeTypes$ALIGNMENT]
getJunctionEdges  <- function(edges) edges[, edgeType != edgeTypes$ALIGNMENT]  

# matchable junctions are compared across reads to accumulate SV evidence
# unmatchable junctions failed one or more quality checks (flanks, bandwidth, adapters, ...) and are ignored
getMatchableJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction & 
        mapQ >= env$MIN_MAPQ & # require confident flanking alignments...
        alnSize >= env$MIN_ALIGNMENT_SIZE & # ... of sufficiently large size
        gapCompressedIdentity >= env$MIN_ALIGNMENT_IDENTITY & # ... and quality, as judged by match to the reference genome
        passedBandwidth == TRUE &
        hasAdapter5 == FALSE &
        hasAdapter3 == FALSE
    ]
}
# preserve junction matchable flag to stratify the quality of the junctions that end up in junction clusters
setMatchableFlag <- function(edges){ 
    isJunction  <- getJunctionEdges(edges)
    isMatchable <- getMatchableJunctions(edges)
    edges[isJunction, matchable := isMatchable[isJunction]] # edge matchable == NA, junction == (TRUE|FALSE)
    edges
}

# fusable junctions are maintained in a single segment, thereby calling a _recurring_ SV
# unfusable junctions split reads into segments and can only be called as singletons
getFusableJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction &  # these filters reject any read not suitable for assembling genome SV paths from multiple reads
        hasAdapter5 == FALSE & # adapters identify nanopore/basecaller artifacts
        hasAdapter3 == FALSE & # each segment side of a read-through junction is kept at this stage, since duplex foldback inversions were collapsed upstream
        (
            passedBandwidth == FALSE |  # we don't consider these to be true junctions, they don't break segments and are ignored during segment matching
            (
                mapQ >= env$MIN_MAPQ & # low-quality alignments can't be trusted to assemble paths
                alnSize >= env$MIN_ALIGNMENT_SIZE & 
                gapCompressedIdentity >= env$MIN_ALIGNMENT_IDENTITY & 
                (
                    nCanonical    > 1 |  # junctions not validated by at least 2 reads on the same strand are likely to be ligation artifacts
                    nNonCanonical > 1    # the required two reads might have been in the same channel (but usually aren't)
                )
            )
        )
    ]
}

#-------------------------------------------------------------------------------------
# segment filters
#-------------------------------------------------------------------------------------

# matchable segments have clustered junctions corresponding to two or more edges (usually in different reads)
# they are thus segments useful for creating SV-driven allelic assemblies
getMatchableSegments <- function(segments){
    segments[, 
        nClusteredJxns    > 0 &  # this segment has at least one clustered junction (really, all of them must)
        nMatchingSegments > 0    # at least one other segment also had at least one of this segment's junctions
    ]
}
