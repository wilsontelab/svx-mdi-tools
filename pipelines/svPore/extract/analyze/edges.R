#-------------------------------------------------------------------------------------
# edges.R carries functions that excute initial parsing and characterization of edges
#-------------------------------------------------------------------------------------

# scan each read to assign properties to individual edges
getBlockN <- function(N, edgeTypes){ # fuse blocks across
    if(N == 3){
        if(edgeTypes[2] %in% splitTypes) 1:3 else 1
    } else {
        blocks <- rle(sapply(edgeTypes, "%in%", splitTypes))$lengths
        unlist(sapply(1:length(blocks), function(i) rep(i, blocks[i]), simplify = FALSE))
    }    
}
checkJunctionFlanks <- function(N, edges){
    passed <- rep(FALSE, N)
    ai <- seq(1, N, 2)
    passed[ai] <- edges[ai, # where alignments are scored to assign their usefulness for scoring high-confidence SVs
        mapQ >= env$MIN_MAPQ & 
        eventSize >= env$MIN_ALIGNMENT_SIZE &
        gapCompressedIdentity >= env$MIN_ALIGNMENT_IDENTITY
    ]
    ji <- seq(2, N, 2)
    passed[ji] <- sapply(ji, function(i) passed[i - 1] && passed[i + 1])
    passed    
}
parseEdgeMetadata <- function(edges){
    message("expanding edge metadata")
    edges <- cbind(
        edges, 
        edges[, parseSignedWindow(node1, 1)], 
        edges[, parseSignedWindow(node2, 2)]
    )
    setkey(edges, qName)

    message("numbering edges and alignment blocks")
    edges[, ":="(
        blockN = getBlockN(.N, edgeType), # number block of alignments on single chrom-strands
        edgeN = 1:.N, # number the edges in each molecule; alignments are odd, junctions are even
        nEdges = .N,    
        cigar = ifelse(is.na(cigar), cigar2, cigar),
        passedFlankCheck = checkJunctionFlanks(.N, .SD)
    ), by = .(qName)]

    # clean up extraneous columns; always use this syntax to drop data.table columns!
    edges[, cigar2 := NULL] 
    edges[, blastIdentity2 := NULL]
    edges[, gapCompressedIdentity2 := NULL]

    ######### restrict tmp file size while developing (or, make this deletion permanent?)
    edges[, cigar := NULL]

    setkey(edges, qName, blockN, edgeN)
    edges
}

# scan alignment blocks to assess their bandwidth
fillBandwidthNodePairs <- function(N){
    i <- seq(2, N, 2)
    x <- as.data.table(expand.grid(leftmost = i, rightmost = i))
    x <- x[rightmost >= leftmost]
    x[, nNodes := rightmost - leftmost + 1]
    x[order(-nNodes)]
} 
checkBlockBandwidth <- function(edges, i1, i2){
    nQryBases <- edges[i2 + 1, xStart] - edges[i1 - 1, xEnd]
    nRefBases <- edges[i2, xEnd] - edges[i1, xStart]
    abs(nQryBases - nRefBases) >= env$MIN_SV_SIZE
}
checkIndelBandwidth <- function(edges){
    message("checking indel bandwith to identify low quality alignment blocks")
    nodePairs <- list()
    edges[, passedBandwidth := {
        if(.N == 1) TRUE # "A" blocks or "T|V" junctions 
        else if(.N == 3) c(TRUE, checkBlockBandwidth(.SD, 2, 2), TRUE) # single "D|U|I" junctions
        else {
            passed <- rep(TRUE, .N) # "ADAIADA" and other complex block paths
            k <- as.character(.N)
            if(is.null(nodePairs[[k]])) nodePairs[[k]] <<- fillBandwidthNodePairs(.N)
            for(j in 1:nrow(nodePairs[[k]])){
                i1 <- nodePairs[[k]][j]$leftmost
                i2 <- nodePairs[[k]][j]$rightmost
                failed <- !checkBlockBandwidth(.SD, i1, i2)
                if(failed) passed[i1:i2] <- FALSE # don't write passes to prevent overwriting prior failures over wider junction spans
            }
            passed
        }
    }, by = .(qName, blockN)]
    edges
}
