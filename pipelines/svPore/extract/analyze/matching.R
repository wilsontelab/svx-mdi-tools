#-------------------------------------------------------------------------------------
# matching.R has functions that support:
#   re-orientation of alignments and junctions to allow proper matching between opposite strands
#   fuzzy SV path matching, i.e., with tolerance for matching nodes in closely adjacent windows due to boundary effects
#-------------------------------------------------------------------------------------

# rectify an edge, i.e., a pair of nodes, to a consistently ordered string representation suitable for sorting and grouping
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0 # always returns the same ordered nodes independent of strand (it is unimportant what that order is)
}
isCanonicalStrandPath <- function(nodePath){
    isCanonicalStrand(c(nodePath[1], nodePath[length(nodePath)]))
}
getCanonicalNodes <- function(edges){
    isAlignment <- getAlignmentEdges(edges)
    x <- edges[, 
        if(isAlignment[.I]) .(
            isCanonical = NA, # alignment canonical is determined per segment, later
            cNode1 = NA_integer_, 
            cNode2 = NA_integer_,
            cChromIndex1 = chromIndex1,
            cChromIndex2 = chromIndex1
        ) else if(isCanonicalStrand(c(node1, node2))) .( 
            isCanonical = TRUE,
            cNode1 = node1, # "c" prefix means "canonical"
            cNode2 = node2,
            cChromIndex1 = chromIndex1,
            cChromIndex2 = chromIndex2
        ) else .(
            isCanonical = FALSE,
            cNode1 = -node2, # reverse the node pair to the opposite strand
            cNode2 = -node1,
            cChromIndex1 = chromIndex2,
            cChromIndex2 = chromIndex1
        ),
        by = c("qName","blockN","edgeN")
    ]
    x[, ":="(
        junctionKey  = paste(cNode1, cNode2, sep = ":")
    )]
    x[, .SD, .SDcols = c("isCanonical","cNode1","cNode2","cChromIndex1","cChromIndex2","junctionKey")] # suitable for cbind'ing to edges
}
# sortCanonical <- function(nodes){
#     isCanonical <- isCanonicalStrand(c(nodes[1], nodes[length(nodes)]))
#     if(!isCanonical) nodes <- -rev(nodes)
#     nodes
# }
# getPathSignature <- function(nodes, reverse = FALSE){
#     isCanonical <- isCanonicalStrand(c(nodes[1], nodes[length(nodes)]))
#     if(!isCanonical) nodes <- -rev(nodes)
#     if(reverse) nodes <- -rev(nodes)
#     paste0(":", paste(nodes, collapse = ":"), ":")
# }

# operators for performing fuzzy and exact node path matching
# handle node pairs or longer node paths
pathTolerance <- 1 # i.e., "fuzzy" allows a one-window adjacency tolerance
pathOffsets <- -pathTolerance:pathTolerance
`%is_fuzzy_match%` <- function(path1, path2){ # core function that checks if two equal-length node sequences are an adjacency match
    all(apply(sapply(pathOffsets, "+", path2) - path1, 1, function(v) min(abs(v))) == 0)
}
`%is_fuzzy_rev_match%` <- function(path1, path2){ # check if two fuzzy paths are reciprocal with respect to strand
    path1 %is_fuzzy_match% -rev(path2)
}
`%is_exact_match%` <- function(path1, path2){ # related function applied after assignment of reference junctions, where exact matches are expected
    identical(path1, path2)
}
`%is_exact_rev_match%` <- function(path1, path2){ # check if two exact paths are exactly reciprocal with respect to strand
    path1 %is_exact_match% -rev(path2)
}
