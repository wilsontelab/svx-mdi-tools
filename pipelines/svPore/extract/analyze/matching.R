#-------------------------------------------------------------------------------------
# matching.R has functions that support re-orientation of alignments and junctions 
# to allow proper matching between opposite strands
#-------------------------------------------------------------------------------------

# rectify an edge, i.e., a pair of nodes, to a consistently ordered string representation suitable for sorting and grouping
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0 # always returns the same ordered nodes independent of strand (it is ~unimportant what that order is)
}
isCanonicalStrandPath <- function(nodePath){
    isCanonicalStrand(c(nodePath[1], nodePath[length(nodePath)]))
}
setCanonicalNodes <- function(edges){
    isAlignment <- getAlignmentEdges(edges)
    x <- edges[, 
        if(isAlignment[.I]) .(
            isCanonical = NA, # alignment canonical is determined per segment, later
            cChromIndex1 = chromIndex1,
            cChromIndex2 = chromIndex1,
            cStrand1 = NA_integer_,
            cStrand2 = NA_integer_,
            cRefPos1 = NA_integer_,
            cRefPos2 = NA_integer_
        ) else if(isCanonicalStrand(c(node1, node2))) .( 
            isCanonical = TRUE,
            cChromIndex1 = chromIndex1, # "c" prefix means "canonical"
            cChromIndex2 = chromIndex2,
            cStrand1 = as.integer(sign(node1)),
            cStrand2 = as.integer(sign(node2)),
            cRefPos1 = xStart,
            cRefPos2 = xEnd 
        ) else .(
            isCanonical = FALSE,
            cChromIndex1 = chromIndex2, # reverse the node pair to the opposite strand
            cChromIndex2 = chromIndex1,
            cStrand1 = as.integer(-sign(node2)),
            cStrand2 = as.integer(-sign(node1)),
            cRefPos1 = xEnd,
            cRefPos2 = xStart 
        ),
        by = c("qName","blockN","edgeN")
    ]
    edges <- merge(edges, x, by = c("qName","blockN","edgeN"), all.x = TRUE)
    setkey(edges, qName, blockN, edgeN)
    isJunction <- getJunctionEdges(edges)
    edges[isJunction, ":="(
        junctionKey  = paste(cChromIndex1, cChromIndex2, cStrand1, cStrand2, cRefPos1, cRefPos2, insertSize, sep = ":")
    )]
    edges
}

# assemble readPaths for reads and segments
getIndexedReadPath <- function(chromSizes, edges, usableJunctions){ # does NOT include outer alignment endpoints
    edges[usableJunctions, {
        iRefPos1 <- if(is.na(iRefPos1)) xStart else iRefPos1 # fallback for real but unmatchable junctions that weren't scored by findMatchingJunctions
        iRefPos2 <- if(is.na(iRefPos2)) xEnd   else iRefPos2
        .(
            node1      = chromSizes[chromIndex1, sign(node1) * (nBasesBefore + iRefPos1 - 1)], # complete description of a single indexed junction 
            insertSize = if(is.na(iInsertSize)) insertSize else iInsertSize, 
            node2      = chromSizes[chromIndex2, sign(node2) * (nBasesBefore + iRefPos2 - 1)]  
        )
    }, by = .(edgeId)][, .SD, .SDcols = c("node1","insertSize","node2")] # complete description of a series of indexed junctions, as a data.table with integer64
}
getCanonicalReadPath <- function(iReadPath, isCanonical){
    if(isCanonical) return(iReadPath)
    iReadPath[, ":="(
        node1 = -node1,
        node2 = -node2
    )]
    rev(iReadPath)[order(.N:1)]
}
getCanonicalPathKey <- function(cOuterNode1, icReadPath, cOuterNode2){
    paste(
        cOuterNode1,
        paste(
            apply(data.table(
                n1 = as.character(icReadPath[[1]]), # must do as.character by column before apply since integer64
                is = as.character(icReadPath[[2]]),
                n2 = as.character(icReadPath[[3]])
            ), 1, paste, collapse = ":"), 
            collapse = "::"
        ),
        cOuterNode2,
        sep = ":::"
    )
}


# # operators for performing node path matching
# isOuterEndpointMatch <- function(cOuterNode1_1, cOuterNode2_1, cOuterNode1_2, cOuterNode2_2){
#     cOuterNode1_1 == cOuterNode1_2 &&
#     cOuterNode2_1 == cOuterNode2_2
# }
# isReadPathMatch <- function(icReadPath_1, icReadPath_2){ # icReadPaths are data.tables of indexed, canonical junction sequences
#     nrow(icReadPath1) == nrow(icReadPath2) && 
#     all(icReadPath1 == icReadPath2)
# }

# pathTolerance <- 1 # i.e., "fuzzy" allows a one-window adjacency tolerance
# pathOffsets <- -pathTolerance:pathTolerance
# `%is_fuzzy_match%` <- function(path1, path2){ # core function that checks if two equal-length node sequences are an adjacency match
#     all(apply(sapply(pathOffsets, "+", path2) - path1, 1, function(v) min(abs(v))) == 0)
# }
# `%is_exact_match%` <- function(path1, path2){ # related function applied after assignment of reference junctions, where exact matches are expected
#     identical(path1, path2)
# }
# `%is_fuzzy_rev_match%` <- function(path1, path2){ # check if two paths are fuzzy reciprocal with respect to strand
#     path1 %is_fuzzy_match% -rev(path2)
# }
# `%is_exact_rev_match%` <- function(path1, path2){ # check if two paths are exactly reciprocal with respect to strand
#     path1 %is_exact_match% -rev(path2)
# }
