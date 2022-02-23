#=====================================================================================
# functions for parsing nodes into networks, i.e., contigs
#-------------------------------------------------------------------------------------
# NB: at present, code does NOT build full networks, it:
#   characterizes individual junctions, i.e., one at a time
#   merges gaps and outer clips into SV sets when they match sequenced junctions
#   records which other junctions are carried in the same molecules as the index junctions
#-------------------------------------------------------------------------------------
# thus:
#   all sequenced junctions with a unique node pair will generate an output SV call
#   gap junctions are either (but not both):
#       present in an SV call keyed by a sequenced junction or a different gap junction, or
#       create an output SV call as the index junction
#   outer clips are either:
#       present in an SV call keyed by a different junction
#       omitted from the output (outer clips alone cannot call a junction)
# where junction calls carry the information needed to follow a network, i.e., junction chain, i.e., allele
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
        abs(diff(pos)) <= env$PURGE_DISTANCE # allows for small alignment differences at junction
    }
})
getMatchingNodes <- function(nodeName, nodeN, nodeClass){
    
    # get nodes on similar chrom-strand within SV distance of query node
    nodeData <- unpackNodeNames(nodeName)
    nodePos  <- nodeData[[3]]
    nodes    <- getNodes('proximity', nodeName, unpackNodeNames = TRUE)
    if(is.null(nodes)) return(NULL)
    nodes[, nodeN := ..nodeN]

    # be more precise about exactly which nodes match the query SV node
    nodeMatches <- if(nodeClass == nodeClasses$GAP){ # query is a gap
        compare <- if(nodeData[[2]] == "L") `>` else `<`
        nodes[, NODE_CLASS == nodeClasses$GAP | # all gaps match based on proximity alone
               compare(pos, nodePos)] # add splits/clips to gaps if they are to the inside
    } else { # query is a split
        compare <- if(nodeData[[2]] == "L") `<` else `>`
        nodes[, (NODE_CLASS == nodeClasses$GAP & compare(pos, nodePos)) | # add gaps to splits if they are to the outside # nolint 
                (NODE_CLASS == nodeClasses$OUTER_CLIP & NODE == nodeName) | # clips must match query split exactly
                (NODE_CLASS == nodeClasses$SPLIT & isSimilarSplit(NODE, nodeName))] # splits allow a fuzzier matching # nolint            
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
    unknownJunctionClass = 4, # rare failure in followSvJunction::junctionMatches
    noUsableNodes = 5,        # error in characterizeSVJunction
    tooFewRefNodes = 6        # error in characterizeSVJunction
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
unknownJunctionClasses <- list()
followSvJunction <- function(edge){

    # parse and track node/junction names
    nodeNames <- c(edge$NODE1, edge$NODE2)
    junctionName <- paste(nodeNames, collapse = ',')
    followedJunctions[junctionName] <<- TRUE

    # get and label nodes consistent with the query junction (e.g., flanking gaps, outer clips)     
    # apply node-level filters
    matchingNodes1 <- getMatchingNodes(edge$NODE1, 1L, edge$nodeClass)
    matchingNodes2 <- getMatchingNodes(edge$NODE2, 2L, edge$nodeClass)
    reject <- applySvNodeFilters(edge, matchingNodes1, matchingNodes2)
    if(!is.null(reject)) return(reject)        
    junctionNodes <- rbind(matchingNodes1, matchingNodes2)    
    junctionNodes[, junctionKey := paste(MOL_ID, JXN_N, sep = ":")]
    
    # further validate that the nodes match on all relevant sides of the query junction
    junctions <- junctionNodes[,
        lapply(.SD, paste, collapse = ","),
        keyby = junctionKey,
        .SDcols = c('NODE', 'NODE_CLASS', 'nodeN')
    ]
    junctionMatches <- mapply(
        function(NODE, NODE_CLASS, nodeN){
            switch(NODE_CLASS,
                '0,0' = nodeN == "1,2",       # GAP, must match on both sides (not both nodes on same side of query) # nolint          
                '1,1' = nodeN == "1,2",       # SPLIT, must match on both sides
                '2'   = TRUE,                 # OUTER_CLIP, one-sided, match to SV already established above
                '0'   = FALSE,                # unpaired nodes (e.g., when one read mapped to a different locus)
                '1'   = FALSE,
                { # failure patterns, record for reporting (future debugging?) and reject the junction
                    if(is.null(unknownJunctionClasses[[NODE_CLASS]])) unknownJunctionClasses[[NODE_CLASS]] <<- 0
                    unknownJunctionClasses[[NODE_CLASS]] <<- unknownJunctionClasses[[NODE_CLASS]] + 1
                    NA 
                }
            )
        },
        junctions[, NODE], junctions[, NODE_CLASS], junctions[, nodeN]
    )    
    if(any(is.na(junctionMatches))) return( rejectJunction('unknownJunctionClass') )
    
    # finalize the list of all nodes that match the query junction
    # must include the query nodes themselves, and potentially other added gaps and outer clips
    # (do not expect new splits to be added, they should have already been merged in previously)
    junctionKeys  <- junctions[junctionMatches, junctionKey]
    junctionNodes <- junctionNodes[junctionKey %in% junctionKeys]    
    
    # alignment MAPQ filter (generally higher than the original group+consensus acceptance filter)
    # applies to two nodes that flank the junction, either a gap or split (but not a clip)
    maxMapQ1 <- junctionNodes[nodeN == 1 & NODE_CLASS != nodeClasses$OUTER_CLIP, max(MAPQ)]
    maxMapQ2 <- junctionNodes[nodeN == 2 & NODE_CLASS != nodeClasses$OUTER_CLIP, max(MAPQ)]
    if(
       !(maxMapQ1 >= env$MIN_MAPQ_BOTH && maxMapQ2 >= env$MIN_MAPQ_BOTH) || # filtered applied to the _best_ gap or split # nolint
       !(maxMapQ1 >= env$MIN_MAPQ_ONE  || maxMapQ2 >= env$MIN_MAPQ_ONE)
    ) return( rejectJunction('failedNodeFilters') )
 
    # mark all of the junctions that match the query as "known" and "followed"
    matchingJunctionNames <- unique(junctions[junctionMatches, NODE])
    matchingJunctionNames <- matchingJunctionNames[matchingJunctionNames != junctionName]
    followedJunctions[matchingJunctionNames] <<- TRUE # thus, non-query junctions added to the evidence list
    knownJunctionNames <- unique(c(junctionName, matchingJunctionNames))

    # get all of the nodes in all of the molecules that matched the query junction
    # includes the query SV nodes
    # adds outermost nodes on those molecules, and maybe additional SV nodes too
    moleculeIds <- junctionNodes[, unique(MOL_ID)]
    moleculeNodes <- getNodes('molecule', moleculeIds, setIsSVJunction = TRUE, unpackNodeNames = TRUE)
    if(is.null(moleculeNodes)) return( rejectJunction('getNodes') )

    # discover whether there are any additional junctions contained in those molecules
    # these are SV junctions that do NOT match the query junction but are also in the query junction molecules
    otherJunctionNames <- if(moleculeNodes[, sum(isSVJunction) > 0]){
        junctionNames <- moleculeNodes[
            isSVJunction == TRUE, 
            list(NODE = paste(NODE, collapse = ",")), # are these guaranteed to be in proper junction-node order?
            keyby = list(MOL_ID, JXN_N)               # depends on data.table sort...
        ][, unique(NODE)] 
        junctionNames[junctionNames %notin% knownJunctionNames]   
    } else character()
     
    # mark the molecule nodes that matched the query junction (to simplify the output to a single nodes list)
    moleculeNodes[, junctionKey := paste(MOL_ID, JXN_N, sep = ":")]
    moleculeNodes[, ':='(IS_JUNCTION_NODE = junctionKey %in% ..junctionKeys,
                         IS_SEED_NODE = NODE %in% ..nodeNames,
                         IS_REPEAT = 0)]
    moleculeNodes <- merge(
        moleculeNodes,
        junctionNodes[, .(junctionKey, ALN_N, NODE_N = nodeN)],
        c('junctionKey', 'ALN_N'),
        all.x = TRUE
    )

    # return our findings
    list(
        rejected = FALSE,
        junctionName = junctionName,
        matchingJunctionNames = matchingJunctionNames, # does not include junctionName
        otherJunctionNames = otherJunctionNames,        
        nodes = moleculeNodes[order(NODE_N)], # "of all the nodes in all the molecules..."
        nJunctionMolecules = length(moleculeIds) # for use during junction reconstruction
    )  
}
#=====================================================================================
