#=====================================================================================
# describe a consensus SV based on sequence data from a set of molecules discovered as evidence
# along the way, purge molecules that are considered false evidence
#-------------------------------------------------------------------------------------
getFractionMatchingBases <- function(seq1, seq2, fixedLength = NULL, side = NULL){
    if(is.null(fixedLength)){
        seq1Len <- nchar(seq1)
        seq2Len <- nchar(seq2)
        if(seq1Len > seq2Len){
            seq1 <- if(side == "L") substr(seq1, 1, seq2Len)
                    else substr(seq1, (seq1Len - seq2Len) + 1, seq1Len)
        } else if(seq2Len > seq1Len) {
            seq2 <- if(side == "L") substr(seq2, 1, seq1Len)
                    else substr(seq2, (seq2Len - seq1Len) + 1, seq2Len)
        }
        fixedLength <- nchar(seq1)
    }
    sum(mapply('==', strsplit(seq1, ''), strsplit(seq2, ''))) / fixedLength
}
#-------------------------------------------------------------------------------------
    # initialize a blank/null junction
#-------------------------------------------------------------------------------------
characterizeSVJunction <- function(nodes){ # nodes pre-sorted by node #
    call <- list(
        rejected = FALSE,
        JXN_TYPE = "",
        JXN_SEQ = NA,   # sequence of the single read or merged reads used to call a junction # nolint
        MERGE_LEN = NA, # number of merged bases, can use as flag for reconstructed junctions
        MICROHOM_LEN = 0,
        MICROHOM_MATCH = NA,
        JXN_BASES = NA,
        refNodes = data.table(), # two nodes that guided construction of the call (whether real or reconstructed)
        nodes = data.table()     # all nodes, with our modifications/additions
    )
   
#-------------------------------------------------------------------------------------
    # examine all SV-flanking molecules for nearly-identical molecule spans
    # aggregate presumed duplicate molecules to prevent them from appearing as falsely independent evidence
#-------------------------------------------------------------------------------------
    junctionKeys <- nodes[, unique(junctionKey)] # MOL_ID:JXN_N   
    nodes[, N_COLLAPSED := 1] # when nothing is purged, every node represents 1 molecule
    if(length(junctionKeys) > 1){ # singleton-molecule SVs need no action

        # generate two identifying outer positions per molecule-junction, one position per each node
        x <- nodes[, {
            list(
                OUT_POS1 = OUT_POS1[1],
                OUT_POS2 = OUT_POS2[1],
                N_SEED_NODE = sum(IS_SEED_NODE), # for determining which molecules to keep
                STRAND_COUNT = STRAND_COUNT1[1] + STRAND_COUNT2[1]
            )
        }, by = "junctionKey"]  

        # calculate the euclidean distance between the endpoints of each pair of molecules
        pairs <- setDT(as.data.frame(t(combn(junctionKeys, 2)), stringsAsFactors = FALSE))
        setnames(pairs, c('junctionKey1', 'junctionKey2'))
        pairs[, dist := apply(pairs, 1, function(ij) x[
            ij,
            dist(.SD[, .(OUT_POS1, OUT_POS2)]),
            on = "junctionKey"
        ])]    

        # for all colliding pairs, keep the molecule with the best read evidence and add other molecules to its counts
        collisionIs <- pairs[, .I[dist < env$PURGE_DISTANCE]]
        if(length(collisionIs) > 0){
            remappedJunctionKeys <- list()
            rejectedMoleculeIds <- integer()
            for(i in collisionIs){     
   
                # determine which molecule of the colliding pair to keep
                bestJ <- which.max(sapply(paste0('junctionKey', 1:2), function(jkCol){
                    x[ # best = highest STRAND_COUNT, but with preference to SEED junctions
                        junctionKey == pairs[i, ..jkCol][[1]],
                        N_SEED_NODE * 1e4 + STRAND_COUNT
                    ]
                }))
                bestJunctionKey <- pairs[i, ..bestJ][[1]]         
                worstJ <- if(bestJ == 1) 2 else 1    
                worstJunctionKey <- pairs[i, ..worstJ][[1]]
  
                # ascend up any prior chain to the best of a set of colliding molecules
                while(!is.null(remappedJunctionKeys[[bestJunctionKey]])) {
                    bestJunctionKey <- remappedJunctionKeys[[bestJunctionKey]]
                }
                remappedJunctionKeys[[worstJunctionKey]] <- bestJunctionKey # the molecule that holds a rejected molecule's counts # nolint

                # add the worst strand counts to the best molecule and count how often we did this
                worstCount1 <- splitGapNodes[junctionKey == worstJunctionKey, STRAND_COUNT1[1]]
                worstCount2 <- splitGapNodes[junctionKey == worstJunctionKey, STRAND_COUNT2[1]]
                nodes[junctionKey == bestJunctionKey, ':='(
                    STRAND_COUNT1 = STRAND_COUNT1 + worstCount1,
                    STRAND_COUNT2 = STRAND_COUNT2 + worstCount2,
                    N_COLLAPSED   = N_COLLAPSED + 1
                )]

                # add to the list of rejected molecules
                rejectedMoleculeIds <- append(rejectedMoleculeIds,
                                              nodes[junctionKey == worstJunctionKey, MOL_ID[1]])
            }   

            # finally, completely purge the presumed duplicate _molecules_ from the nodes list
            nodes <- nodes[!(MOL_ID %in% rejectedMoleculeIds)]
        }
    }    

#-------------------------------------------------------------------------------------
    # select one molecule as a reference for characterizing the junction
    # make sure it has the same signature as the original seed junction for the network
#-------------------------------------------------------------------------------------
    usableNodes <- nodes[, NODE_CLASS == nodeClasses$SPLIT] # splits can always completely sequence a junction, prefer them # nolint
    if(!any(usableNodes)) usableNodes <- TRUE               # gaps can only approximately locate a junction
    # for splits and gaps, prefer seed junctions, i.e., molecules of the type used to build the SV network    
    # for splits, prefer longest clip on the shorter-clip side (promotes long molecules with central junctions)
    # for gaps, prefer molecules with the longest clips on either side
    minMaxFn <- if(nodes[usableNodes, NODE_CLASS][1] == nodeClasses$SPLIT) min else max
    bestJxnKey <- nodes[ 
        usableNodes,
        .(
            N_SEED_NODE   = sum(IS_SEED_NODE),                
            minMaxClipLen = minMaxFn(CLIP_LEN)
        ), 
        by = junctionKey
    ][
      order(-N_SEED_NODE, -minMaxClipLen)
    ][ 
       1,
       junctionKey
    ]
    refNodes <- nodes[usableNodes & junctionKey == bestJxnKey] # always exactly two nodes that help build the call set 
    call$JXN_TYPE <- refNodes[1, JXN_TYPE]
    if(call$JXN_TYPE == junctionTypes$UNKNOWN) return(call)

#-------------------------------------------------------------------------------------
    # ensure that the reference nodes are properly paired and ordered
#-------------------------------------------------------------------------------------
    if(refNodes[, paste0(nodeN, collapse = ",") != "1,2"]) {
        return(list(rejected = TRUE, reason = rejectionReasons$tooFewRefNodes))
    } 

#-------------------------------------------------------------------------------------
    # add any outer clips that match inner clips at the reference nodes (whether reference in a split or a gap)
#-------------------------------------------------------------------------------------
    nodes[, IS_RC := FALSE] # whether or not characterizeSVJunction applied RC to node sequences
    isInversion <- refNodes[1, side] == refNodes[2, side]
    flipNode    <- if(refNodes[1, side] == 'R') 1L else 2L
    for(i in 1:2){
        if(refNodes[i, CLIP_LEN == 0]) next # unclipped gap nodes aren't informative for outer clip matching
        outerClipNodes <- getNodes('outer_clips', refNodes[i, NODE], unpackNodeNames = TRUE)
        if(is.null(outerClipNodes)) next
        outerClipNodes[, ':='(
            junctionKey = paste(MOL_ID, JXN_N, sep = ":"),
            nodeN = i,
            IS_SEED_NODE = FALSE,
            N_COLLAPSED = 1,
            IS_RC = FALSE
        )]

        # flip inversion outer clips as needed to also put them on the canonical strand
        if(isInversion && flipNode == i){        
            outerClipNodes[, ':='(
                SEQ      = rc(SEQ),
                CLIP_SEQ = rc(CLIP_SEQ),
                CIGAR    = rc_cigar(CIGAR),
                IS_RC    = TRUE
            )]
        } 

        # merge inner and outer node evidence
        nodes <- rbind(nodes, outerClipNodes) # TODO: sort clips in by nodeN ??
    }

#-------------------------------------------------------------------------------------
    # attempt to merge gap+clip-only junctions by using clipped nodes from different molecules
    # could be inner or outer clips of any length, but each partner must be clipped to correctly set 'pos'
#-------------------------------------------------------------------------------------
    if(refNodes[1, NODE_CLASS] == nodeClasses$SPLIT) { 
        call$JXN_SEQ <- refNodes[1, SEQ] 
    } else if(nodes[, length(unique(MOL_ID)) > 1]){ # only relevant when there are multiple molecules but no splits
        
        # get the two nodes with the greatest potential for overlap with the other side
        # NB: merge nodes can include outer clips
        mergeNodes <- do.call('rbind', lapply(1:2, function(i){
            d <- nodes[nodeN == i]
            j <- if(refNodes[i, side] == 'L') d[, which.max(pos + CLIP_LEN)]
                                         else d[, which.min(pos - CLIP_LEN)]
            d[j] 
        }))

        # find the best overlap score over all possible merge registers
        if(mergeNodes[, all(CLIP_LEN > 0)] && # both nodes must be clipped, so that node positions characterize a junction # nolint
           mergeNodes[1, MOL_ID] != mergeNodes[2, MOL_ID]){ #  we've already tried to merge nodes from the same molecule # nolint
            maxOverlap <- mergeNodes[, min(nchar(SEQ))]
            if(env$MIN_MERGE_OVERLAP <= maxOverlap){
                end1 <- nchar(mergeNodes[1, SEQ])        
                overlapLengths <- env$MIN_MERGE_OVERLAP:maxOverlap
                mergeScores <- unlist(mclapply(overlapLengths, function(overlapLength){
                    o1 <- substr(mergeNodes[1, SEQ], end1 - overlapLength + 1, end1)
                    o2 <- substr(mergeNodes[2, SEQ], 1, overlapLength)
                    getFractionMatchingBases(o1, o2, overlapLength)
                }, mc.cores = env$N_CPU))
                iAtMax <- which.max(mergeScores) 
    
                # if a proper overlap/merge is found, reconstruct the inferred SV junction
                if(mergeScores[iAtMax] >= env$MIN_MERGE_DENSITY){
                    call$JXN_SEQ <- paste0(
                        substr(mergeNodes[1, SEQ], 1, end1 - overlapLengths[iAtMax]),
                        mergeNodes[2, SEQ]
                    )                    
                    call$MERGE_LEN <- overlapLengths[iAtMax]                
                    refNodes <- mergeNodes
                }                
            }
        }
    }

#-------------------------------------------------------------------------------------
    # purge clip molecules
    #   from unsequenced, gap-only junctions; they cannot be validated as meaningful
    #   that don't match the assembled junction on the unaligned side
#-------------------------------------------------------------------------------------
    clipIs <- nodes[, .I[NODE_CLASS == nodeClasses$OUTER_CLIP]]
    if(length(clipIs) > 0){
        discardClipIs <- if(is.na(call$JXN_SEQ)){
            clipIs # result therefore leaves only gap molecules
        } else {
            na.omit(sapply(clipIs, function(i){
                nodeN <- nodes[i, nodeN]
                x <- getFractionMatchingBases(nodes[i, CLIP_SEQ],
                                              refNodes[nodeN, CLIP_SEQ],
                                              side = refNodes[nodeN, side])
                if(x < env$MIN_MERGE_DENSITY) i else NA
            }))
        }       
        if(length(discardClipIs) > 0) {
            nodes <- nodes[!(MOL_ID %in% nodes[discardClipIs, unique(MOL_ID)])] # remove clips as evidence
        }
    }

#-------------------------------------------------------------------------------------
    # characterize the junction when possible
#-------------------------------------------------------------------------------------
    if(!is.na(call$JXN_SEQ)){
        
        # positions within the complete junction sequence as bounded by refNode1 and refNode2 SEQ limits
        seq1Len <- refNodes[1, nchar(SEQ)]
        readPos1 <- seq1Len - refNodes[1, CLIP_LEN]
        readPos2 <- if(is.na(call$MERGE_LEN)){ # sequenced junction
            refNodes[2, CLIP_LEN] + 1
        } else { # reconstructed junction
            seq1Len - call$MERGE_LEN + refNodes[2, CLIP_LEN] + 1
        }

        # parse microhomology vs. inserted bases
        call$MICROHOM_LEN <- readPos1 - readPos2 + 1
        if(call$MICROHOM_LEN > 0){ # microhomology
            call$JXN_BASES <- substr(call$JXN_SEQ, readPos2, readPos1)
        } else if(call$MICROHOM_LEN < 0){ # inserted bases
            call$JXN_BASES <- substr(call$JXN_SEQ, readPos1 + 1, readPos2 - 1)
        } else { # a blunt joint
            call$JXN_BASES <- ""
        }          
    }
    
#-------------------------------------------------------------------------------------
    # add a flag to nodes to indicate which nodes guided the call construction
#-------------------------------------------------------------------------------------
    nodes[, IS_REF_NODE := FALSE]
    rnKeys <- refNodes[, paste(MOL_ID, JXN_N, ALN_N, sep = ":")]
     nKeys <-    nodes[, paste(MOL_ID, JXN_N, ALN_N, sep = ":")]
    nodes[nKeys %in% rnKeys, IS_REF_NODE := TRUE]
    
#-------------------------------------------------------------------------------------
    # return our result
    # a "call" is a set of molecules supporting a specifically characterized SV junction
    # include refNodes for use in constructing the output table
#-------------------------------------------------------------------------------------
    call$refNodes <- refNodes    
    call$nodes    <- nodes
    call
}
#=====================================================================================
