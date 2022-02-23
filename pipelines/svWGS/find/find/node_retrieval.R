
# functions for gathering information on nodes and edges
# requires data.table

#=====================================================================================
# load required data and indices on script launch
#-------------------------------------------------------------------------------------
# different ways of querying node relationships
#-------------------------------------------------------------------------------------
nodeFiles <- list()
for(by in c('proximity',  # find all nodes within SV distance of a query node (step #1)  
            'molecule')){ # find all nodes clamied by a query molecule (step #2)  

    # make the nodes file ready for indexed reading
    type  <- paste0('nodes_by_', by)
    file  <- paste(type, 'txt', sep = ".")
    file  <- paste(env$SHM_DIR_WRK, file, sep = "/")
    index <- paste(file, 'index', sep = ".")
    nodeFiles[[by]] <- list(
        # connection to the file that was placed en bloc in RAM drive by find_junctions.sh
        con   = file(file, open = "rb"),
        # index to the file, loaded into R
        index = fread(file = index, sep = "\t", header = TRUE, stringsAsFactors = FALSE)   
    )
    
    # standardize the index's key
    setnames(nodeFiles[[by]]$index, 1, 'key')
    if(typeof(nodeFiles[[by]]$index$key) != 'character')
        nodeFiles[[by]]$index[, key := as.character(key)]
    setindex(nodeFiles[[by]]$index, key)
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
    x[[3]] <- as.integer(x[[3]])
    x
}
# load the edges we will analyze as a data.table
edges <- do.call('rbind', lapply(c(
    'sequenced', # SV junctions that were sequenced outright by split alignments (processed first)
    'gap'        # SV junctions in the unsequenced bases between two reads; pos is variable
), function(edgeType){     
    x <- fread(
        file = paste(env$COMPILE_PREFIX, paste0(edgeType, '_junctions'), 'gz', sep = "."),
        sep = "\t",        
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = names(compile$edges),
        colClasses = unname(unlist(compile$edges))
    )
    x[, c('chrom1', 'side1', 'pos1')] <- unpackNodeNames(x$NODE1)
    x[, c('chrom2', 'side2', 'pos2')] <- unpackNodeNames(x$NODE2)    
    edgeClass <- edgeClasses[[edgeType]]
    x[, ':='(edgeType  = ..edgeType,
             nodeClass = ..edgeClass,
             size      = abs(pos2 - pos1))]
    #x[, edgePairing := paste(chrom1, side1, chrom2, side2, sep=":")] # get segFault when tried to parallelize by edgePairing # nolint
    x[order(-COUNT_DISTINCT)] # thus, use most frequent splits first
}))
#=====================================================================================

#=====================================================================================
# random access to sets of nodes by various query types
#-------------------------------------------------------------------------------------
getNodeFileChunks <- Vectorize(function(by, keyValue){
    chunk <- nodeFiles[[by]]$index[as.character(keyValue), on = "key"]
    seek(nodeFiles[[by]]$con, where = chunk$offset, origin = "start", rw = "r")
    readChar(nodeFiles[[by]]$con, chunk$size)
})
getNodes <- function(by, keyValues, unpackNodeNames=FALSE, setIsSVJunction=FALSE){
             
    # use index file to retrieve a data block
    x <- tryCatch({    
        paste0(getNodeFileChunks(by, keyValues), collapse = "")
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
    if(by == 'molecule') x[order(-NODE_CLASS, NODE)] # splits then gaps; data.table radix sort should match linux sort
    else x        
}
#=====================================================================================
