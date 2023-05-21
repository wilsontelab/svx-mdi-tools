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
getAllJunctions <- getJunctionEdges  

# real junctions are compared across reads to find duplex repetition
# they represent all true reference discontinuities, including low-quality and adapter-chimeric junctions,
# but not including low-bandwidth junctions that we ignore as "unreal"
getRealJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction & 
        passedBandwidth == TRUE
    ]
}

# matchable junctions are compared across reads to accumulate SV evidence
# unmatchable junctions failed one or more quality checks (flanks, bandwith, adapters, ...) and are ignored
getMatchableJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction & 
        passedFlankCheck == TRUE & # the junctions we are still willing to call after single-read analysis
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
            passedBandwidth == FALSE |      # we don't consider these to be true junctions, they don't break segments and are ignored during segment matching
            (
                passedFlankCheck == TRUE &  # low-quality alignments can't be trusted to assemble paths
                (
                    nCanonical    > 1 |  # junctions not validated by at least 2 reads on the same strand are likely to be ligation artifacts
                    nNonCanonical > 1    # see notes in duplex.R
                )
            )
        )
    ]
}

# singleton junctions are high quality junctions nominated as potential ultra-rare SVs
# they are analyzed without attempting to assemble genome segment paths around them
# many singleton junctions are artifacts, e.g., due to ligation, but the list can be mined for true SVs
getSingletonJunctions <- function(edges){
    isQualityJunction <- getMatchableJunctions(edges)
    edges[, 
        isQualityJunction & 
        nJunctionInstances == 1 # may have been a duplex detection
    ]
}

#-------------------------------------------------------------------------------------
# segment filters
#-------------------------------------------------------------------------------------

# matchable segments have matchable junctions corresponding to two or more edges (usually in different reads)
# they are thus segments useful for creating SV-driven allelic assemblies
getMatchableSegments <- function(segments){
    segments[, 
        nMatchableJxns > 0 &  # this segment has at least one matchable junction (really, all of them must)
        nMatchingSegments > 0 # at least one other segment also had at least one of this segment's junctions
    ]
}
