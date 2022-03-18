
# functions for gathering information on nodes and edges

#=====================================================================================
# load required data and indices on script launch
#-------------------------------------------------------------------------------------
# nodes, i.e., single breakpoints in the traversal graph
#-------------------------------------------------------------------------------------
nodeFiles <- list()
for(type in c('nodes_by_proximity',  # to find all nodes within SV distance of a query node
              'outer_clips')){       # to find all outer clips matching a split node       

    # make the nodes file ready for indexed reading
    file  <- paste(type, 'txt', sep = ".")
    file  <- paste(env$SHM_DIR_WRK, file, sep = "/")
    index <- paste(file, 'index', sep = ".")
    nodeFiles[[type]] <- list(
        # connection to the file that was placed en bloc in RAM drive by find_junctions.sh
        con   = file(file, open = "rb"),
        # index to the file, loaded into R
        index = fread(file = index, sep = "\t", header = TRUE, stringsAsFactors = FALSE)   
    )
    
    # standardize the index's key, always handled as a string
    setnames(nodeFiles[[type]]$index, 1, 'key')
    if(typeof(nodeFiles[[type]]$index$key) != 'character')
        nodeFiles[[type]]$index[, key := as.character(key)]
    setindex(nodeFiles[[type]]$index, key)
}
#-------------------------------------------------------------------------------------
# two types of two-sided junctions (i.e., SV edges between nodes), in priority order
#-------------------------------------------------------------------------------------
# extract chrom, side, pos from node sort keys, e.g., 12:R:32538
unpackNodeNames <- function(nodeNames){ 
    x <- as.data.frame(
        matrix(unlist(strsplit(nodeNames, '\\W')), ncol = 3, byrow = TRUE),
        stringsAsFactors = FALSE
    )
    x[[3]] <- as.integer(x[[3]]) # TODO: treat chromI as integer?
    x
}
# load the edges we will analyze as a data.table
edges <- fread(
    file = paste(env$COMPILE_PREFIX, 'junction_edges.gz', sep = "."),
    sep = "\t",        
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = names(compile$edges),
    colClasses = unname(unlist(compile$edges))
)
edges[, c('chrom1', 'side1', 'pos1')] <- unpackNodeNames(edges$NODE1)
edges[, c('chrom2', 'side2', 'pos2')] <- unpackNodeNames(edges$NODE2) 
edges[, size := abs(pos2 - pos1)]
edges <- edges[order(-NODE_CLASS, -COUNT_DISTINCT)] # thus, use most frequent splits first
#=====================================================================================

#=====================================================================================
# random access to sets of nodes by various query types
#-------------------------------------------------------------------------------------
getNodeFileChunks <- Vectorize(function(type, keyValue){
    chunk <- nodeFiles[[type]]$index[as.character(keyValue), on = "key"]
    seek(nodeFiles[[type]]$con, where = chunk$offset, origin = "start", rw = "r")
    readChar(nodeFiles[[type]]$con, chunk$size)
})
getNodes <- function(type, keyValues, unpackNodeNames = FALSE, setIsSVJunction = FALSE){
             
    # use index file to retrieve a data block
    x <- tryCatch({    
        paste0(getNodeFileChunks(type, keyValues), collapse = "")
    }, error = function(e){ # typically because non-canonical chrom present in query but not target nodes
        NULL
    })
    if(is.null(x)) return(NULL)

    # build a data table of all unique matching nodes
    # concatenate the text block for each of a series of query nodes, then read as table
    x <- fread(
        text = x,
        sep = "\t",        
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = names(compile$nodes),
        colClasses = unname(unlist(compile$nodes))
    )

    # parse/sort node information as needed
    if(unpackNodeNames) x[, c('chrom', 'side', 'pos')] <- unpackNodeNames(x$NODE)
    if(setIsSVJunction) x[,
        isSVJunction := (NODE_CLASS == nodeClasses$SPLIT) |
                        (NODE_CLASS == nodeClasses$GAP & JXN_TYPE != junctionTypes$PROPER)
    ]
    x        
}
#=====================================================================================
