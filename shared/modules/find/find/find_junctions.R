
# find and characterize molecule groups that call source junctions

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(abind)
library(parallel)
library(data.table)
#-------------------------------------------------------------------------------------
# load and parse environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(rUtilDir, 'workflow.R')
loadEnvVars(list(
    string = c(
        'ACTION_DIR',
        'GENOME',
        'COMPILE_PREFIX',
        'FIND_PREFIX', 
        'SHM_DIR_WRK'
    ),
    integer = c(
        'N_CPU',
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
    ),
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
#options(warn=2)
#-------------------------------------------------------------------------------------
# source R scripts
source(rUtilDir, 'utilities.R')
#-------------------------------------------------------------------------------------
sequenceUtilDir <- file.path(rUtilDir, 'sequence')
source(file.path(sequenceUtilDir, 'general.R'))
source(file.path(sequenceUtilDir, 'IUPAC.R'))
source(file.path(sequenceUtilDir, 'smith_waterman.R'))
#-------------------------------------------------------------------------------------
genomeUtilDir <- file.path(rUtilDir, 'genome')
source(file.path(genomeUtilDir, 'general.R'))
source(file.path(genomeUtilDir, 'faidx.R'))
#-------------------------------------------------------------------------------------
findStepDir <- file.path(env$ACTION_DIR, 'find')
sourceScript(findStepDir, 'column_definitions.R')
sourceScript(findStepDir, 'node_retrieval.R')
sourceScript(findStepDir, 'network.R')
sourceScript(findStepDir, 'analyze_junctions.R')
#-------------------------------------------------------------------------------------
# initialize genome
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
processEdge <- function(i){
    
    # parse the query junction
    junction <- edges[i, ]
    junctionName <- paste(junction[, NODE1], junction[, NODE2], sep = ',')
    if(!is.null(followedJunctions[[junctionName]])) return( followedJxn ) # edge already handled with a previous junction # nolint

    # build a network of molecules containing the junction
    network <- followSvJunction(junction)
    if(network$rejected) return( rejectEdge(network$reason) )
    
    # characterize the assembled junction
    call <- characterizeSVJunction(network$nodes, network$nJunctionMolecules)
    if(call$rejected) return( rejectEdge(call$reason) )

    # get attributes of molecules and junction nodes in the network
    molAgg_first <- call$nodes[ 
        IS_JUNCTION_NODE == TRUE,
        lapply(.SD, function(v) v[1]),
        by = MOL_ID,
        .SDcols = molAggCols_first
    ]
    molAgg_collapse <- call$nodes[ 
        IS_JUNCTION_NODE == TRUE,
        lapply(.SD, paste, collapse = ":"),
        by = MOL_ID,
        .SDcols = molAggCols_collapse
    ]
    
    # assign an SV unique identifier
    SV_ID <<- SV_ID + 1
    call$nodes[, SV_ID := ..SV_ID]    

    # append all nodes in all matching molecules to disk
    # store indices for fast node retrieval
    index <- printNodes(SV_ID, call$nodes)

    # parse the filterable table of all SV junctions
    paste(  
        SV_ID,
        network$junctionName,
        paste(network$matchingJunctionNames, collapse = "::"),
        paste(network$otherJunctionNames,    collapse = "::"),
        #-------------
        junction$TARGET_CLASS, # NB: _not_ the same as the _molecule_ target class attached to nodes
        call$JXN_TYPE,        
        #-------------
        call$refNodes[1, chrom], # TODO: convert chrom indices back to string chroms
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
        paste(nchar(molAgg_first$SEQ), collapse = ","),
        paste(molAgg_collapse$UMI,  collapse = ","),
        paste(molAgg_collapse$MAPQ, collapse = ","),
        #-------------
        index$CHUNK_OFFSET,
        index$CHUNK_SIZE,
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
previouslyFollowed <- svTable[, is.na(TARGET_CLASS)]
rejectedNodes <- function(reason) svTable[, TARGET_CLASS == as.character(reason)]
rejectedGetNodes       <- rejectedNodes(rejectionReasons$getNodes)
rejectedMatchingNodes  <- rejectedNodes(rejectionReasons$matchingNodes)
rejectedNodeFilter     <- rejectedNodes(rejectionReasons$failedNodeFilters)
rejectedUnknownClass   <- rejectedNodes(rejectionReasons$unknownJunctionClass)
rejectedNoUsableNodes  <- rejectedNodes(rejectionReasons$noUsableNodes)
rejectedTooFewRefNodes <- rejectedNodes(rejectionReasons$tooFewRefNodes)
svTable <- svTable[!previouslyFollowed &                   
                   !rejectedGetNodes &
                   !rejectedMatchingNodes &
                   !rejectedNodeFilter &    
                   !rejectedUnknownClass &
                   !rejectedNoUsableNodes &
                   !rejectedTooFewRefNodes, ]
save(svTable, file = paste(env$FIND_PREFIX, 'structural_variants', 'RData', sep = "."))
#-------------------------------------------------------------------------------------
# report some more stats
#-------------------------------------------------------------------------------------
nRejected <- function(v) sum(v, na.rm = TRUE)
reportStat(nrow(svTable),                     "SV junctions reported")
reportStat(nRejected(rejectedGetNodes),       "edges rejected because getNodes threw an error")
reportStat(nRejected(rejectedMatchingNodes),  "edges rejected because inversion nodes matched each other")
reportStat(nRejected(rejectedNodeFilter),     "edges rejected due to node filter options (e.g., MAPQ)")
reportStat(nRejected(rejectedUnknownClass),   "edges rejected due to unknown junction class in followSvJunction")
reportStat(nRejected(rejectedNoUsableNodes),  "edges rejected with no usable nodes in characterizeSVJunction")
reportStat(nRejected(rejectedTooFewRefNodes), "edges rejected with !=2 ref nodes in characterizeSVJunction")
#=====================================================================================
