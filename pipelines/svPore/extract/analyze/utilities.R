#-------------------------------------------------------------------------------------
# svPore analyzeNodes support functions
#-------------------------------------------------------------------------------------
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0
}
getPathSignature <- function(NODES){
    isCanonical <- isCanonicalStrand(c(NODES[1], NODES[length(NODES)]))
    if(!isCanonical) NODES <- -rev(NODES)
    paste(NODES, collapse = ":")
}
parseSignedWindow <- function(window, side) {
    strand <- ifelse(window > 0, "+", "-")
    window <- abs(window)
    chromIndex <- bitwShiftR(window, 24)
    dt <- data.table(
        chromIndex  = chromIndex,
        # chrom       = chromIndex ? $revChromIndex{$chromIndex} : "?",
        chrom       = unlist(revChromIndex[chromIndex]),
        windowIndex = bitwAnd(window, 2**24 - 1),
        strand      = strand
    )
    names(dt) <- paste0(names(dt), side)
    dt
}
