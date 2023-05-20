#-------------------------------------------------------------------------------------
# svPore analyze general support functions
#-------------------------------------------------------------------------------------

# load input data
loadEdges <- function(type) {
    message(paste("loading", type, "edges"))
    fread(
        switch(
            type,
            sv  = env$EDGES_SV_FILE,    # file with reads that had split alignments
            tmp = env$EDGES_TMP_FILE    # file with ~10K reads that did not have split alignments
        ),
        col.names   = edgesCols,
        colClasses  = edgesColClasses,
        sep = "\t",
        quote = ""
    )
}
loadReads <- function(){
    message("loading read sequences")
    fread(
        env$SEQUENCES_FILE, # file with source base-level data for reads in env$EDGES_SV_FILE (J lines) + env$EDGES_TMP_FILE (A lines)
        col.names   = readsCols,
        colClasses  = readsColClasses,
        sep = "\t",
        quote = ""
    )
}

# nodes are codified into an integer for streamlined comparison
# this function expands integer nodes out to chrom/strand/pos
parseSignedWindow <- function(chromWindows, nodes, side) {
    strands <- ifelse(nodes > 0, "+", "-")
    genomeIndices <- abs(nodes)
    dt <- cbind(
        chromWindows[genomeIndices + 1, .(chrom, chromIndex, windowIndex)],        
        data.table(
            genomeIndex = genomeIndices,             
            strand      = strands         
        )
    )        
    names(dt) <- paste0(names(dt), side)
    dt
}
