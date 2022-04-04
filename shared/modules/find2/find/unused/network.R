#=====================================================================================
# functions for parsing SV node relationships
#-------------------------------------------------------------------------------------
# NB: at present, code does NOT build full networks, i.e., allelic contigs, it:
#   characterizes individual junctions, i.e., one at a time
#   merges gaps and outer clips into SV sets when they match sequenced junctions
#-------------------------------------------------------------------------------------
# thus:
#   all sequenced junctions with a unique node pair will generate an output SV call
#   gap junctions are either (but not both):
#       present in an SV call keyed by a sequenced junction or a different gap junction, or
#       create an output SV call as the index junction
#   outer clips are either:
#       present in an SV call keyed by a different junction
#       omitted from the output (outer clips alone cannot call a junction)
#=====================================================================================

#=====================================================================================
# get all individual alignment nodes that match a specific query node with respect to chrom/side/pos
# i.e., that are consistent with tracking the same SV junction on one side
#-------------------------------------------------------------------------------------
isSimilarSplit <- Vectorize(function(node1, node2){
    # NB: nodes provided to this function were similar by proximity, so chrom and side must match
    # purpose is to check splits more aggressively by position than done previously 
    node1 == node2 ||
    {
        pos <- unpackNodeNames(c(node1, node2))[[3]]
        abs(diff(pos)) <= env$PURGE_DISTANCE # allow for small alignment differences at sequenced junctions between molecules # nolint
    }
})
getMatchingNodes <- function(nodeName, otherNodeName, nodeN, nodeClass){

    # get nodes on similar chrom-strand within SV distance of query node
    nodeData      <- unpackNodeNames(nodeName)
    otherNodeData <- unpackNodeNames(otherNodeName)
    indexName <- paste(nodeName, otherNodeData[[1]], otherNodeData[[2]], sep=":") # node:partner
    nodes <- getNodes('nodes_by_proximity', indexName)
    if(is.null(nodes)) return(NULL)
    nodes[, NODE_N := ..nodeN]

    # be more precise about exactly which nodes match the query SV node
    nodeMatches <- if(nodeClass == nodeClasses$GAP){ # query is a gap
        compare <- if(nodeData[[2]] == "L") `>` else `<`
        nodes[, NODE_CLASS == nodeClasses$GAP | # all gaps match based on proximity alone
                compare(pos, nodeData[[3]])]    # add splits to gaps if they are to the inside
    } else { # query is a split
        compare <- if(nodeData[[2]] == "L") `<` else `>`
        nodes[, (NODE_CLASS == nodeClasses$GAP & compare(pos, nodeData[[3]])) | # add gaps to splits if they are to the outside # nolint 
                isSimilarSplit(NODE, nodeName)]                                 # splits allow fuzzy matching # nolint            
    }
    nodes[nodeMatches, ]
}
#=====================================================================================

#=====================================================================================
# build a network of all nodes (SV and not) and edges (junctions and alignments)
# that fully describe all source molecules that contain a single specific SV junction
# the resulting network may not be fully closed if a _different_ SV also exists in those molecules
#-------------------------------------------------------------------------------------
# establish rejection criteria and flags, listed in the order they are applied
rejectionReasons <- list(
    getNodes = 1,             # rare unexpected failures to extract the nodes from the indexed files (why does this happen?) # nolint
    matchingNodes = 2,        # the two nodes in a junction are the same as each other (rare small inversions in bad regions) # nolint
    failedNodeFilters = 3,    # user didn't request processing of the junction (e.g., insufficent MAPQ, etc.),
    tooFewRefNodes = 4        # error in characterizeSVJunction
)
rejectJunction <- function(reason) list(rejected = TRUE, reason = rejectionReasons[[reason]])
applySvNodeFilters <- function(junction, nodes1, nodes2){ # filters that apply to individual alignment nodes

    # check on the existence of data for the nodes
    if(
       is.null(nodes1) ||
       is.null(nodes2)
    ) return(list(rejected = TRUE, reason = rejectionReasons$getNodes))

    # ensure that the nodes are not the same, i.e., call a true SV
    if(
        nodes1[, any(NODE == junction$NODE2)] ||
        nodes2[, any(NODE == junction$NODE1)]
    ) return(list(rejected = TRUE, reason = rejectionReasons$matchingNodes))
    return(NULL) # NULL = success, i.e., all filters passed
}
#-------------------------------------------------------------------------------------
followedJunctions <- list() # mark junctions refused or added to a previous SV call set to prevent them from creating a new call set # nolint
followSvJunction <- function(edge){

    # parse and track node/junction names
    nodeNames <- c(edge$NODE1, edge$NODE2)
    junctionName <- paste(nodeNames, collapse = ',')
    followedJunctions[[junctionName]] <<- TRUE

    # get and label all junction nodes consistent with the query junction   
    matchingNodes1 <- getMatchingNodes(edge$NODE1, edge$NODE2, 1L, edge$NODE_CLASS)
    matchingNodes2 <- getMatchingNodes(edge$NODE2, edge$NODE1, 2L, edge$NODE_CLASS)

    # apply node-level rejection filters    
    reject <- applySvNodeFilters(edge, matchingNodes1, matchingNodes2)
    if(!is.null(reject)) return(reject)

    # restrict to source molecules with node on both sides of the junction
    # note that outer clip molecules have not yet been added
    jxnMolIds <- intersect(matchingNodes1$MOL_ID, matchingNodes2$MOL_ID) # always includes the query junction
    # TODO: here and/or later could add a rejection criteria based on junction molecule coverage
    matchingNodes1 <- matchingNodes1[MOL_ID %in% jxnMolIds]
    matchingNodes2 <- matchingNodes2[MOL_ID %in% jxnMolIds]
    
    # alignment MAPQ filter (generally higher than the original group+consensus acceptance filter)
    # applies to two nodes that flank the junction, either a gap or split (but not a clip)
    maxMapQ1 <- matchingNodes1[, max(MAPQ)]
    maxMapQ2 <- matchingNodes2[, max(MAPQ)]
    if(
       !(maxMapQ1 >= env$MIN_MAPQ_BOTH && maxMapQ2 >= env$MIN_MAPQ_BOTH) || # filtered applied to the _best_ gap or split # nolint
       !(maxMapQ1 >= env$MIN_MAPQ_ONE  || maxMapQ2 >= env$MIN_MAPQ_ONE)
    ) return( rejectJunction('failedNodeFilters') )

    # assemble the list of all junctions matching the query and mark them as followed
    junctionNodes <- rbind(matchingNodes1, matchingNodes2)
    junctionNodes[, junctionKey := paste(MOL_ID, JXN_N, sep = ":")]
    junctions <- junctionNodes[,
        lapply(.SD, paste, collapse = ","),
        keyby = junctionKey,
        .SDcols = c('NODE')
    ]
    matchingJunctionNames <- unique(junctions[, NODE])
    followedJunctions[matchingJunctionNames] <<- TRUE # thus, non-query junctions added to the evidence list

    # mark the nodes that matched the query junction
    junctionNodes[, IS_SEED_NODE := 0 + (NODE %in% ..nodeNames)]

    # return our findings
    list(
        rejected = FALSE,
        junctionName = junctionName, # the query junction
        matchingJunctionNames = matchingJunctionNames, # includes junctionName  
        nodes = junctionNodes
    )  
}
#=====================================================================================
