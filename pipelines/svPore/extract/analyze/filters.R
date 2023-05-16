#-------------------------------------------------------------------------------------
# filters.R returns filters of edges and segments based on their assessed properties
# all functions return a logical vector of the same length as the provided data.table
#-------------------------------------------------------------------------------------

# all alignment or junction edges
getAlignmentEdges <- function(edges) edges[, edgeType == edgeTypes$ALIGNMENT]
getJunctionEdges  <- function(edges) edges[, edgeType != edgeTypes$ALIGNMENT]  

# matchable junctions are compared across reads to accumulate SV evidence
# unmatchable junctions failed one or more quality checks (flanks, bandwith, adapters, ...) and are ignored
getMatchableJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction & 
        passedFlankCheck == TRUE & # in other words, the junctions we are still willing to call after single-read analysis
        passedBandwidth == TRUE &
        hasAdapter == FALSE
    ]
}

# fusable junctions are maintained in a single segment, thereby calling a _recurring_ SV
# unfusable junctions split reads into segments
getFusableJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction &  # these filters reject any read not suitable for assembling genome SV paths from multiple reads
        !hasAdapter & # adapters identify nanopore/basecaller artifacts
        (
            passedBandwidth == FALSE | # we don't consider these to be true junctions, they don't breaks segments and are ignored during segment matching
            (
                passedFlankCheck == TRUE & # low-quality alignments can't be trusted to assemble paths
                (
                    nCanonical    > 1 | # two detections on one strand guarantee that a junction was observed in two independent DNA molecules
                    nNonCanonical > 1
                )
            )
        )
    ]
}

# matchable segments have matchable junctions corresponding to two or more edges (usually in different reads)
# they are thus segments useful for chaining SV junctions
getMatchableSegments <- function(segments){
    segments[, 
        nMatchableJxns > 0 &  # this segment has at least one matchable junction (really, all of them must)
        nMatchingSegments > 0 # at least one other segment also had at least one of this segment's junctions
    ]
}

