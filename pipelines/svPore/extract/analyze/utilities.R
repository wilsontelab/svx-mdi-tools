#-------------------------------------------------------------------------------------
# svPore analyze general support functions
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
