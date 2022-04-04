# functions for gathering information on nodes and edges

#=====================================================================================
# load required data and indices on script launch
#-------------------------------------------------------------------------------------
# nodes, i.e., single breakpoints in the traversal graph
#-------------------------------------------------------------------------------------
nodeFiles <- list()
connectNodeFile <- function(sample, type){

    # make the nodes file ready for indexed reading
    file  <- paste(type, 'txt', sep = ".")
    file  <- paste(env$SHM_DIR_WRK, sample, file, sep = "/")
    index <- paste(file, 'index', sep = ".")
    x <- list(
        # connection to the file that was placed en bloc in RAM drive by find_junctions.sh
        con   = file(file, open = "rb"),
        # index to the file, loaded into R
        index = fread(file = index, sep = "\t", header = TRUE, stringsAsFactors = FALSE)   
    )
    
    # standardize the index's key, always handled as a string
    setnames(x$index, 1, 'key')
    if(typeof(x$index$key) != 'character')
        x$index[, key := as.character(key)]
    setindex(x$index, key)
    x
}
#-------------------------------------------------------------------------------------
# extract chrom, side, pos from node sort keys, e.g., 12:R:32538
unpackNodeNames <- function(nodeNames){ 
    x <- as.data.frame(
        matrix(unlist(strsplit(nodeNames, '\\W')), ncol = 3, byrow = TRUE),
        stringsAsFactors = FALSE
    )
    x[[1]] <- as.integer(x[[1]])
    x[[3]] <- as.integer(x[[3]])
    x
}
#=====================================================================================

#=====================================================================================
# random access to sets of nodes by various query types
#-------------------------------------------------------------------------------------
getNodeFileChunks <- Vectorize(function(sample, type, keyValue){
    if(is.null(nodeFiles[[sample]][[type]])) nodeFiles[[sample]][[type]] <<- connectNodeFile(sample, type)
    chunk <- nodeFiles[[sample]][[type]]$index[as.character(keyValue), on = "key"]
    if(is.na(chunk$offset)) return("")
    seek(nodeFiles[[sample]][[type]]$con, where = chunk$offset, origin = "start", rw = "r")
    readChar(nodeFiles[[sample]][[type]]$con, chunk$size)
})
getNodes <- function(sample, type, keyValues){
    # concatenate the text block for each of a series of query nodes, then read as table
    x <- paste0(getNodeFileChunks(sample, type, keyValues), collapse = "")
    if(nchar(x) == 0) return(NULL)
    fread(
        text = x,
        sep = "\t",        
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = names(compile$nodes),
        colClasses = unname(unlist(compile$nodes))
    )
}
getOuterClips <- function(sample, nodes, svIdx){
    x <- getNodes(sample, 'outer_clips', nodes)
    if(is.null(x)) return(NULL)
    x[, ':='( # convert nodes to molecules with duplicated node
        SAMPLE           = sample,
        NODE_2           = NODE_1,
        CLIP_LEN_2       = CLIP_LEN_1,
        CLIP_SEQ_2       = CLIP_SEQ_1,
        #---------------
        FLAG_2           = FLAG_1,
        POS_2            = POS_1,
        MAPQ_2           = MAPQ_1,
        CIGAR_2          = CIGAR_1,
        SEQ_2            = SEQ_1,
        ALN_N_2          = ALN_N_1,        
        #---------------
        UMI_2            = UMI_1,
        #---------------
        groupIndex       = 0,
        AMBIGUOUS        = 0,
        N_COLLAPSED      = 1,
        jxnName          = paste(NODE_1, NODE_1, sep = ","), # the junction edge called by a molecule (could be >1 per molecule)
        jxnKey           = paste(SAMPLES[sample], MOL_ID, JXN_N, sep = ":") # ID for the source junction edge (each with 2 nodes)
    )]
    x[, c('chrom1', 'side1', 'pos1')] <- unpackNodeNames(x$NODE_1)
    x[, c('chrom2', 'side2', 'pos2')] <- x[, c('chrom1', 'side1', 'pos1')]
    x[, ':='( # convert nodes to molecules with duplicated node
        svIndex          = svIdx,
        sampleSvIndex    = paste(SAMPLES[sample], svIdx, sep=":")
    )]
    x        
}
#=====================================================================================
