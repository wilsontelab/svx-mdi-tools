# matching.R has functions that support:
#   fuzzy SV path matching, i.e., with tolerance for matching nodes in closely adjacent windows
#   partial SV path matching, i.e., that reports on molecule overlap/nesting as well as duplex/identity

# rectify a junction edge, i.e., a pair of nodes, to a consistently ordered string representation suitable for sorting and grouping
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0
}
getCanonicalNodes <- function(edges){
    x <- edges[, 
        if(isCanonicalStrand(c(node1, node2))) .( # junctions (but not full paths) can be re-ordered to a canonical strand for matching
            isCanonical = TRUE,
            cNode1 = node1, # "c" prefix means "canonical"
            cNode2 = node2,
            cChromIndex1 = chromIndex1,
            cChromIndex2 = chromIndex2
        ) else .(
            isCanonical = FALSE,
            cNode1 = -node2,
            cNode2 = -node1,
            cChromIndex1 = chromIndex2,
            cChromIndex2 = chromIndex1
        ),
        by = c("segmentName","blockN","edgeN")
    ]
    x[, ":="(
        junctionKey  = paste(cNode1, cNode2, sep = ":")
    )]
    x[, .SD, .SDcols = c("isCanonical","cNode1","cNode2","cChromIndex1","cChromIndex2","junctionKey")]
}

# operators for performing fuzzy and exact junction matching
pathTolerance <- 1
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
