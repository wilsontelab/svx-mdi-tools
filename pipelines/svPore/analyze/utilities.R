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
# loadReads <- function(){
#     message("loading read sequences")
#     fread(
#         # env$SEQUENCES_FILE, # file with source base-level data for reads in env$EDGES_SV_FILE + env$EDGES_TMP_FILE
#         cmd = paste("cut -f 1-2", env$SEQUENCES_FILE),
#         # col.names   = readsCols,
#         # colClasses  = readsColClasses,
#         col.names   = c("qName","qSeq"),
#         colClasses  = c("character","character"),
#         sep = "\t",
#         quote = ""
#     )[, .SD, .SDcols = c("qName","qSeq")]
# }

# nodes are codified into an integer64 for streamlined comparison
# this function expands integer nodes out to chrom/strand/pos
parseSignedNodes <- function(chromSizes, nodes, side) {
    genomeIs <- abs(nodes)
    chromIs <- Vectorize(function(i) which(chromSizes$nBasesBefore >= i)[1])(genomeIs) # sapply does not work with integer64!
    dt <- chromSizes[chromIs, .(
        chrom = chrom, 
        chromIndex = chromIndex,
        refPos = genomeIs - nBasesBefore,
        strand = ifelse(nodes > 0, "+", "-")
    )]    
    names(dt) <- paste0(names(dt), side)
    dt
}
