#-------------------------------------------------------------------------------------
# svPore analyze support functions
#-------------------------------------------------------------------------------------
# load input data
loadEdges <- function(type) {
    message(paste("loading", type, "edges"))
    fread(
        if(type == "sv") env$EDGES_SV_FILE else env$EDGES_TMP_FILE,
        col.names = edgesCols,
        colClasses = edgesColClasses,
        sep = "\t",
        quote = ""
    )
}
loadReads <- function(){
    message("loading read sequences")
    fread(
        env$SEQUENCES_FILE,
        col.names = readsCols,
        colClasses = readsColClasses,
        sep = "\t",
        quote = ""
    )
}

# expand integer nodes out to chrom/strand/pos
parseSignedWindow <- function(window, side) {
    strand <- ifelse(window > 0, "+", "-")
    window <- abs(window)
    chromIndex <- bitwShiftR(window, 24)
    dt <- data.table(
        chromIndex  = chromIndex,
        chrom       = unlist(revChromIndex[chromIndex]),
        windowIndex = bitwAnd(window, 2**24 - 1),
        strand      = strand
    )
    names(dt) <- paste0(names(dt), side)
    dt
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


# junction edges are pairs of nodes flanking the junction itself
# used as first quick check for potential segment matches, i.e., if !any(getJunctionMatch), getSegmentMatch is not needed
# junctionMatchTyes <- list(
#     NONE            = 0,
#     SAME_STRAND     = 1,
#     OPPOSITE_STRAND = 2
# )
# getJunctionMatch <- function(edge1, edge2){ 
#     if(edge1 %is_rev_path% edge2)  junctionMatchTyes$OPPOSITE_STRAND 
#     else if(edge1 %is_path% edge2) junctionMatchTyes$SAME_STRAND
#     else junctionMatchTyes$NONE
# }
# getJunctionMatches <- function(segmentName1, blockN1, edgeN1, edgeType_, chromPair_, edge1){ # check a query edge against all other edges, including itself as a sanity check
#     x <- junctionsToMatch[
#         edgeType == edgeType_ & # single junctions must have a common edgeType (whereas molecule mights have different composite pathTypes)
#         chromPair == chromPair_ & # restrict the search space to the chromosome pair for speed
#         segmentName != segmentName1 # don't match molecules to themselves        
#     ]
#     if(nrow(x) == 0) return(NULL)
#     x[,
#         .( matchType = getJunctionMatch(edge1, c(node1, node2)) ), 
#         by = .(segmentName, blockN, edgeN)
#     ][
#         matchType > junctionMatchTyes$NONE, 
#         .(
#             segmentName = segmentName1,
#             blockN = blockN1,
#             edgeN = edgeN1,
#             matchingSegments = list(segmentName)
#         )
#     ]
# }

