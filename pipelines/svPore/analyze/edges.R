#-------------------------------------------------------------------------------------
# edges.R carries functions that excute initial parsing and characterization of edges
#-------------------------------------------------------------------------------------

# # scan each read to assign properties to individual edges
# getBlockN <- function(N, edgeTypes){ # fuse alignments into blocks across junctions on the same chrom+strand
#     if(N == 3){
#         if(edgeTypes[2] %in% splitTypes) 1:3 else 1
#     } else {
#         blocks <- rle(sapply(edgeTypes, "%in%", splitTypes))$lengths
#         unlist(sapply(1:length(blocks), function(i) rep(i, blocks[i]), simplify = FALSE))
#     }    
# }
# checkJunctionFlanks <- function(N, edges){
#     passed <- rep(FALSE, N)
#     ai <- seq(1, N, 2)

#     # score alignments to assess their usefulness for scoring high-confidence SVs
#     passed[ai] <- edges[ai, 
#         mapQ >= env$MIN_MAPQ & # require confident alignments...
#         eventSize >= env$MIN_ALIGNMENT_SIZE & # ... of sufficiently large size
#         gapCompressedIdentity >= env$MIN_ALIGNMENT_IDENTITY # ... and quality, as judged by match to the reference genome
#     ]

#     # pass junctions if _both_ of their flanking junctions passed
#     ji <- seq(2, N, 2)
#     passed[ji] <- sapply(ji, function(i) passed[i - 1] && passed[i + 1])
#     passed    
# }
# parseEdgeMetadata <- function(edges, chromSizes){
#     message("expanding edge metadata")

#     # edges <- cbind(
#     #     edges, 
#     #     edges[, parseSignedNodes(chromSizes, node1, 1)], 
#     #     edges[, parseSignedNodes(chromSizes, node2, 2)]
#     # )
#     setkey(edges, qName)

#     message("numbering edges and alignment blocks")
#     edges[, edgeId := 1:.N] # number edges uniquely across all reads for junction matching
#     edges[, ":="(
#         blockN = getBlockN(.N, edgeType), # number block of alignments on single chrom-strands
#         edgeN = 1:.N, # number the edges in each molecule; alignments are odd, junctions are even
#         nEdges = .N,    
#         passedFlankCheck = checkJunctionFlanks(.N, .SD)
#     ), by = .(qName)]

#     setkey(edges, qName, blockN, edgeN)
#     edges
# }

# scan alignment blocks to assess their bandwidth
# lower quality base stretches in nanopore reads can lead to large indels where size(I) ~= size(D)
# that are almost never of interest as true SVs
fillBandwidthNodePairs <- function(N){
    i <- seq(2, N, 2)
    x <- as.data.table(expand.grid(leftmost = i, rightmost = i))
    x <- x[rightmost >= leftmost]
    x[, nNodes := rightmost - leftmost + 1]
    x[order(-nNodes)]
} 
checkBlockBandwidth <- function(edges, i1, i2){
    nQryBases <- edges[i2 + 1, qStart] - edges[i1 - 1, qEnd]
    nRefBases <- abs(edges[i2 + 1, node1] - edges[i1 - 1, node2]) #  edges[i2, rEnd] - edges[i1, rStart]
    abs(nQryBases - nRefBases) >= env$MIN_SV_SIZE
}
checkIndelBandwidth <- function(edges){
    message("checking block bandwiths to identify low quality base stretches with false indels")
    nodePairs <- list()
    edges[edgeType != edgeTypes$SPACER, passedBandwidth := {
        if(.N == 1) TRUE # "A" blocks or "T|V" junctions 
        else if(.N == 3) c(TRUE, checkBlockBandwidth(.SD, 2, 2), TRUE) # single "D|U|I" junctions
        else {
            passed <- rep(TRUE, .N) # "ADAIADA" and other complex block paths
            k <- as.character(.N)
            if(is.null(nodePairs[[k]])) nodePairs[[k]] <<- fillBandwidthNodePairs(.N)
            for(j in 1:nrow(nodePairs[[k]])){
                i1 <- nodePairs[[k]][j]$leftmost
                i2 <- nodePairs[[k]][j]$rightmost
                if(!checkBlockBandwidth(.SD, i1, i2)) passed[i1:i2] <- FALSE # don't write passes to prevent overwriting prior failures over wider junction spans
            }
            passed
        }
    }, by = .(edgeSetN, blockN)]
    edges
}
