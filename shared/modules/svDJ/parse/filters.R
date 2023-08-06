#-------------------------------------------------------------------------------------
# filters.R returns filters of edges and segments based on their assessed properties
# all functions return a logical vector of the same length as the provided data.table
#-------------------------------------------------------------------------------------
# this svDJ file extends svPore/filters.R
#-------------------------------------------------------------------------------------

# fusable junctions are maintained in a single segment representing a single, independent amplicon instance
getNonChimericJunctions <- function(edges){
    isJunction <- getJunctionEdges(edges)
    edges[, 
        isJunction &
        isChimeric == FALSE
    ]
}

# V(D)J junction grouping is simplified by only considering segments with one junction
# even allowing one _matchable_ junction creates problems due to low quality reads
getSingletonMatchableJunctions <- function(edges){
    isMatchable <- getMatchableJunctions(edges)
    edges[, 
        isMatchable &  
        nTotalJunctions == 1 # not nKeptJunctions, which can pass low quality reads
    ]
}
