#-------------------------------------------------------------------------------------
# attempt to merge gap+clip-only junctions by using clipped nodes from different molecules
# could be inner or outer clips of any length, but each partner must be clipped to correctly set 'pos'
#-------------------------------------------------------------------------------------
mergeGapJunctions <- function(svIdx){

    # process only relevant when there are multiple independent molecules but no splits
    jxnMols <- jxnMols[svIndex == svIdx] 
    if(
        jxnMols[, any(NODE_CLASS == nodeClasses$SPLIT)] || 
        jxnMols[, length(unique(MOL_ID)) == 1]
    ){
        return(jxnMols)
    }   

    # get the two nodes with the greatest potential for overlap with the other side
    # NB: merge nodes can include outer clips
    refMol <- jxnMols[IS_REFERENCE == 1]
    mergeMols <- do.call('rbind', lapply(1:2, function(i){
        j <- if(i == 1){
            if(refMol[, side1] == 'L') 
                jxnMols[, which.max(ifelse(NODE_1 == "*", NA, pos1 + CLIP_LEN_1))]
            else 
                jxnMols[, which.min(ifelse(NODE_1 == "*", NA, pos1 - CLIP_LEN_1))]
        } else {
            if(refMol[, side2] == 'L') 
                jxnMols[, which.max(ifelse(NODE_2 == "*", NA, pos2 + CLIP_LEN_2))]
            else 
                jxnMols[, which.min(ifelse(NODE_2 == "*", NA, pos2 - CLIP_LEN_2))]
        }
        jxnMols[j] 
    }))
    if(nrow(mergeMols) != 2) return(jxnMols)
    mergeNodes <- rbind(
        mergeMols[1, .(
            MOL_ID = MOL_ID,
            CLIP_LEN = CLIP_LEN_1,
            SEQ = SEQ_1
        )],
        mergeMols[2, .(
            MOL_ID = MOL_ID,
            CLIP_LEN = CLIP_LEN_2,
            SEQ = SEQ_2
        )]
    )

    # find the best overlap score over all possible merge registers
    if(mergeNodes[, all(CLIP_LEN > 0)] && # both nodes must be clipped, so that node positions characterize a junction # nolint
        mergeNodes[1, MOL_ID] != mergeNodes[2, MOL_ID]){ #  we've already tried to merge nodes from the same molecule # nolint
        maxOverlap <- mergeNodes[, min(nchar(SEQ))]
        if(env$MIN_MERGE_OVERLAP <= maxOverlap){
            end1 <- nchar(mergeNodes[1, SEQ])        
            overlapLengths <- env$MIN_MERGE_OVERLAP:maxOverlap
            mergeScores <- unlist(lapply(overlapLengths, function(overlapLength){
                o1 <- substr(mergeNodes[1, SEQ], end1 - overlapLength + 1, end1)
                o2 <- substr(mergeNodes[2, SEQ], 1, overlapLength)
                getFractionMatchingBases(o1, o2, overlapLength)
            }))
            iAtMax <- which.max(mergeScores) 

            # if a proper overlap/merge is found, reconstruct the inferred SV junction
            if(mergeScores[iAtMax] >= env$MIN_MERGE_DENSITY){
                jxnMols[IS_REFERENCE == 1, ':='(
                    JXN_SEQ = paste0(
                        substr(mergeNodes[1, SEQ], 1, end1 - overlapLengths[iAtMax]),
                        mergeNodes[2, SEQ]
                    ),
                    MERGE_LEN = overlapLengths[iAtMax],
                    RECONSTRUCTED = 1
                )]               
            }                
        }
    }
    jxnMols
}
#-------------------------------------------------------------------------------------
# describe a consensus SV based on sequence data from a set of molecules discovered as evidence
#-------------------------------------------------------------------------------------
characterizeSvJunction <- function(seq1, seq2, fixedLength = NULL, side = NULL){
    # initialize a blank/null junction    
    # call <- list(
    #     rejected = FALSE,
    #     JXN_TYPE = "",
    #     JXN_SEQ = NA,   # sequence of the single read or merged reads used to call a junction # nolint
    #     MERGE_LEN = NA, # number of merged bases, can use as flag for reconstructed junctions
    #     MICROHOM_LEN = 0,
    #     MICROHOM_MATCH = NA,
    #     JXN_BASES = NA,
    #     refNodes = data.table(), # two nodes that guided construction of the call (whether real or reconstructed)
    #     nodes = data.table()     # all nodes, with our modifications/additions
    # )
    # # call$JXN_TYPE <- refMol[1, JXN_TYPE]
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
    # add a flag to indicate which _reference_ nodes guided the call construction
    # this is distinct from the _seed_ nodes used to initialize the SV network
#-------------------------------------------------------------------------------------
    nodes[, IS_REF_NODE := 0]
    rnKeys <- refNodes[, paste(MOL_ID, JXN_N, ALN_N, sep = ":")]
     nKeys <-    nodes[, paste(MOL_ID, JXN_N, ALN_N, sep = ":")]
    nodes[nKeys %in% rnKeys, IS_REF_NODE := 1]
    
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
