#-------------------------------------------------------------------------------------
# edges.R carries functions that excute initial parsing and characterization of edges
#-------------------------------------------------------------------------------------

# scan reads to assess their bandwidth
# lower quality base stretches in nanopore reads can lead to large indels where size(I) ~= size(D)
# or to false, low-quality alignments of those bases to ectopic locations
# these things are almost never of interest as true SVs
fillBandwidthNodePairs <- function(N){
    i <- seq(2, N, 2)
    x <- as.data.table(expand.grid(leftmost = i, rightmost = i))
    x <- x[rightmost >= leftmost]
    x[, nNodes := rightmost - leftmost + 1]
    x[order(-nNodes)]
} 
checkNodePairBandwidth <- function(edges, i1, i2){
    nQryBases <- edges[i2 + 1, qStart] - edges[i1 - 1, qEnd]
    nRefBases <- abs(edges[i2 + 1, node1] - edges[i1 - 1, node2]) # thus, bandwidth filter succeeds for T and V
    abs(nQryBases - nRefBases) >= env$MIN_SV_SIZE
}
checkReadBandwidth <- function(edges){
    message("checking read bandwidths to identify low quality stretches with false indels or alignments")
    nodePairs <- list()
    edges[, passedBandwidth := {
        if(.N == 3) c(NA, checkNodePairBandwidth(.SD, 2, 2), NA) # single junctions
        else {
            passed <- rep(NA, .N) # multi-junction paths
            k <- as.character(.N)
            if(is.null(nodePairs[[k]])) nodePairs[[k]] <<- fillBandwidthNodePairs(.N)
            for(j in 1:nrow(nodePairs[[k]])){
                i1 <- nodePairs[[k]][j]$leftmost
                i2 <- nodePairs[[k]][j]$rightmost
                if(!checkNodePairBandwidth(.SD, i1, i2)) passed[i1:i2] <- FALSE # don't write passes to prevent overwriting prior failures over wider junction spans
            }
            passed
        }
    }, by = .(readI)]
    edges
}
