
# find and characterize molecule groups that call source junctions

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(parallel)
library(data.table)
library(jsonlite)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'DATA_NAME',
        'ACTION_DIR',
        'GENOME',
        'COMPILE_PREFIX',
        'FIND_PREFIX', 
        'SHM_DIR_WRK',
        'GENOME_CHROMS'
    ),
    integer = c(
        'N_CPU',
        'READ_LEN',
        'MAX_TLEN',
        'MIN_MAPQ_ONE',
        'MIN_MAPQ_BOTH',
        'MIN_SV_SIZE',
        'SV_SIZE_FACTOR',
        'PURGE_DISTANCE',
        'MIN_MERGE_OVERLAP',
        'ON_TARGET'
    ),
    double = c(
        'MIN_MERGE_DENSITY'
    )
))
writeEnvJson(env$FIND_PREFIX) # for use by compare and/or Shiny app
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general', 'IUPAC', 'smith_waterman'))
sourceScripts(file.path(rUtilDir, 'genome'),     c('general', 'chroms', 'faidx'))
sourceScripts(file.path(env$ACTION_DIR, 'find'), c('column_definitions', 'node_retrieval', 'network', 'analyze_junctions'))
#-------------------------------------------------------------------------------------
# initialize genome
setCanonicalChroms()
loadFaidx(env$SHM_DIR_WRK)
faidx_padding <- round(env$MAX_TLEN * 1.5, 0) # sufficient to contain any source molecule span
#=====================================================================================

#=====================================================================================
# print indexed file of evidence nodes
#-------------------------------------------------------------------------------------
nodesOffset <- 0
allNodesFile <- paste(env$FIND_PREFIX, 'all_nodes', 'txt', sep = ".")
allNodesCols <- names(find$all_nodes)
printNodes <- function(svId, nodes){

    # print the file with data on all nodes in the molecules that matched an SV
    append <- nodesOffset > 0 # thus creates/overwrites the file on first written SV
    x <- paste(apply(nodes[, ..allNodesCols], 1, paste, collapse = "\t"), collapse = "\n")
    x <- paste0(x, "\n")
    cat(x, file = allNodesFile, append = append)
    
    # return indexing information
    offset <- nodesOffset
    size   <- nchar(x)
    nodesOffset <<- nodesOffset + size       
    list(CHUNK_OFFSET = offset, CHUNK_SIZE = size)
}
#=====================================================================================

#=====================================================================================
# assemble extracted SV molecules into SV calls potentially including many supporting molecules
#-------------------------------------------------------------------------------------
# molecule-level values that are the same on both nodes of a junction
molAggCols_first <- c( 
    'NODE_CLASS',
    'IS_MERGED',    
    'IS_DUPLEX',
    'STRAND_COUNT1',
    'STRAND_COUNT2',
    'SHARED_PROPER',
    'SEQ' # only valid if IS_MERGED = TRUE (otherwise, TLEN is NA to us at this point)
)
# alignment-level values, i.e., potentially different on the two nodes of a junction
molAggCols_collapse <- c( 
    'UMI',
    'MAPQ'
)
#-------------------------------------------------------------------------------------
# apply junction-level filters
#-------------------------------------------------------------------------------------
message("finding SV junctions as sets of anomalous edges between alignment nodes")
nAllEdges <- nrow(edges)

# on target filter
if(env$ON_TARGET != 0){ # no target filter
    edges <- switch(env$ON_TARGET,
        edges[TARGET_CLASS != '--'],     # at least one end in padded target
        edges[!grepl('-', TARGET_CLASS)] # both ends in padded target
    )    
}
nOnTargetEdges  <- nrow(edges)
nOffTargetEdges <- nAllEdges - nOnTargetEdges

# SV size filter
if(env$SV_SIZE_FACTOR > 0){
    minSize <- env$SV_SIZE_FACTOR * env$MAX_TLEN
    edges <- edges[JXN_TYPE == "T" | size >= minSize]
} else if(env$MIN_SV_SIZE > 0){
    edges <- edges[JXN_TYPE == "T" | size >= env$MIN_SV_SIZE]
}
nAnalyzed <- nrow(edges)
nTooSmall <- nOnTargetEdges - nAnalyzed

# report some stats
reportStat(nAllEdges,       "input edges")
reportStat(nOffTargetEdges, "edges rejected by on-target requirements")
reportStat(nTooSmall,       "edges rejected because the SV is too small")
reportStat(nAnalyzed,       "edges subjected to further processing")
#-------------------------------------------------------------------------------------
# process candidate junctions
#-------------------------------------------------------------------------------------
nSvTableCol <- length(find$structural_variants)
rejectEdge  <- function(reason) paste(rep(reason, nSvTableCol), collapse = "\t")
followedJxn <- rejectEdge(NA)
SV_ID <- 0
messageIncrement <- 1000
messageThreshold <- 0
processEdge <- function(i){
    
    # parse the query junction
    junction <- edges[i, ]
    junctionName <- paste(junction[, NODE1], junction[, NODE2], sep = ',')
    if(!is.null(followedJunctions[[junctionName]])) return( followedJxn ) # edge already handled with a previous junction # nolint

    # provide log file feedback
    nFollowed <- length(followedJunctions)
    if(nFollowed >= messageThreshold) {
        message(paste(
            paste0("working on edge ", i),
            paste0("(", nFollowed, " edges already aggregated into SVs)")
        ))
        messageThreshold <<- floor(nFollowed / messageIncrement) * messageIncrement + messageIncrement
    }

    # build a network of molecules containing the junction
    network <- followSvJunction(junction)
    if(network$rejected) return( rejectEdge(network$reason) )
    
    # characterize the assembled junction
    call <- characterizeSVJunction(network$nodes)
    if(call$rejected) return( rejectEdge(call$reason) )

    # get attributes of molecules and junction nodes in the network
    molAgg_first <- call$nodes[,
        lapply(.SD, function(v) v[1]),
        by = MOL_ID,
        .SDcols = molAggCols_first
    ]
    molAgg_collapse <- call$nodes[,
        lapply(.SD, paste, collapse = ":"),
        by = MOL_ID,
        .SDcols = molAggCols_collapse
    ]
    
    # assign an SV unique identifier
    SV_ID <<- SV_ID + 1
    call$nodes[, SV_ID := ..SV_ID]   

    # append all nodes in all matching molecules to disk
    # store indices for fast node retrieval
    call$nodes[, IS_REPEAT := 0]
    index <- printNodes(SV_ID, call$nodes)

    # revert chrom indices back to string names
    call$refNodes[, chrom := unlist(revChromIndex[chrom])]

    # parse the filterable table of all SV junctions
    paste(  
        SV_ID,
        network$junctionName,
        paste(network$matchingJunctionNames, collapse = "::"),
        #-------------
        junction$TARGET_CLASS, # NB: _not_ the same as the _molecule_ target class attached to nodes
        call$JXN_TYPE,        
        #-------------
        call$refNodes[1, chrom],
        call$refNodes[1, side],  
        call$refNodes[1, pos],
        call$refNodes[2, chrom],
        call$refNodes[2, side],  
        call$refNodes[2, pos],
        #-------------
        call$JXN_SEQ,
        call$MERGE_LEN,
        faidx_padding, # sample-specific (not junction-specific, but need value for junction display)
        getRefSeq_padded(call$refNodes[1, chrom], call$refNodes[1, pos], faidx_padding), # never rc'ed yet
        getRefSeq_padded(call$refNodes[2, chrom], call$refNodes[2, pos], faidx_padding),       
        #-------------
        call$MICROHOM_LEN,
        call$MICROHOM_MATCH,
        call$JXN_BASES,        
        if(call$JXN_TYPE == "T") 0 else abs(call$refNodes[2, pos] - call$refNodes[1, pos]),          
        #-------------
        molAgg_first[, .N],
        molAgg_first[, sum(NODE_CLASS == nodeClasses$SPLIT)],
        molAgg_first[, sum(NODE_CLASS == nodeClasses$GAP)],
        molAgg_first[, sum(NODE_CLASS == nodeClasses$OUTER_CLIP)],
        molAgg_first[NODE_CLASS != nodeClasses$OUTER_CLIP, sum(IS_DUPLEX)],
        molAgg_first[                                    , sum(IS_DUPLEX)],
        molAgg_first[NODE_CLASS != nodeClasses$OUTER_CLIP, sum(STRAND_COUNT1 + STRAND_COUNT2)],
        molAgg_first[                                    , sum(STRAND_COUNT1 + STRAND_COUNT2)],
        molAgg_first[NODE_CLASS != nodeClasses$OUTER_CLIP, sum(SHARED_PROPER)],
        molAgg_first[                                    , sum(SHARED_PROPER)], 
        #-------------
        paste(molAgg_first$IS_MERGED,  collapse = ","), # for distribution plots mostly
        paste(nchar(molAgg_first$SEQ), collapse = ","),
        paste(molAgg_collapse$UMI,     collapse = ","),
        paste(molAgg_collapse$MAPQ,    collapse = ","),
        #-------------
        index$CHUNK_OFFSET, # for rapid indexed retrieval of nodes supporting this SV
        index$CHUNK_SIZE,
        #-------------
        env$DATA_NAME, # begin to prepare for subsequent inter-sample comparison
        #-------------
        sep = "\t"
    )  
}
#-------------------------------------------------------------------------------------
# use fread to help parse and assemble the final output table (as RData file for use by Shiny)
svTable <- fread(
    text = paste0( sapply(1:nAnalyzed, processEdge), collapse = "\n" ),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE, 
    col.names  = names(find$structural_variants),
    colClasses = unname(unlist(find$structural_variants))
)
#-------------------------------------------------------------------------------------
# finalize and write the output SV summary table
#-------------------------------------------------------------------------------------
message("writing SV summary table")
previouslyFollowed <- svTable[, is.na(TARGET_CLASS)] # TARGET_CLASS unimportant, all fields have same value when rejected
rejectedNodes <- function(reason) svTable[, TARGET_CLASS == as.character(reason)]
rejectedGetNodes       <- rejectedNodes(rejectionReasons$getNodes)
rejectedMatchingNodes  <- rejectedNodes(rejectionReasons$matchingNodes)
rejectedNodeFilter     <- rejectedNodes(rejectionReasons$failedNodeFilters)
rejectedTooFewRefNodes <- rejectedNodes(rejectionReasons$tooFewRefNodes)
svTable <- svTable[!previouslyFollowed &                   
                   !rejectedGetNodes &
                   !rejectedMatchingNodes &
                   !rejectedNodeFilter &    
                   !rejectedTooFewRefNodes, ]
# save(svTable, file = paste(env$FIND_PREFIX, 'structural_variants', 'RData', sep = "."))
outFile <- paste(env$FIND_PREFIX, 'structural_variants', 'gz', sep = ".")
outFile <- gzfile(outFile, "w")
write.table(
    svTable, 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = FALSE
)
close(outFile)
#-------------------------------------------------------------------------------------
# report some more stats
#-------------------------------------------------------------------------------------
nRejected <- function(v) sum(v, na.rm = TRUE)
reportStat(nrow(svTable),                     "SV junctions reported")
reportStat(nRejected(rejectedGetNodes),       "edges rejected because getNodes threw an error")
reportStat(nRejected(rejectedMatchingNodes),  "edges rejected because inversion nodes matched each other")
reportStat(nRejected(rejectedNodeFilter),     "edges rejected due to node filter options (e.g., MAPQ)")
reportStat(nRejected(rejectedTooFewRefNodes), "edges rejected with !=2 ref nodes in characterizeSVJunction")
#=====================================================================================
