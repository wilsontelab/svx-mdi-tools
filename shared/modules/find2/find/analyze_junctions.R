#=====================================================================================
# describe a consensus SV based on sequence data from a set of molecules discovered as evidence
# along the way, purge molecules that are considered false evidence
#-------------------------------------------------------------------------------------
# assess whether an unaligned sequences, e.g., a clip, matches expectations of aligned reads
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
characterizeSVJunction <- function(svIdx){ # junction nodes (no clips) pre-sorted by NODE_N
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
    # call$JXN_TYPE <- refMol[1, JXN_TYPE]


    # select one molecule as a reference for characterizing the junction
    #   prefer splits, which can always completely sequence a junction  
    #   for splits, prefer longest clip on the shorter-clip side (promotes long molecules with central junctions)
    #   for gaps, prefer molecules with the longest clips on either side
    jxnMols <- jxnMols[svIndex == svIdx]    
    nodes <- rbind( # some actions act per-node, split them apart here
        jxnMols[, .(jxnKey = jxnKey, NODE_N = 1, NODE_CLASS = NODE_CLASS, CLIP_LEN = CLIP_LEN_1)],
        jxnMols[, .(jxnKey = jxnKey, NODE_N = 2, NODE_CLASS = NODE_CLASS, CLIP_LEN = CLIP_LEN_2)]
    )
    usable <- nodes[, NODE_CLASS == nodeClasses$SPLIT]
    if(!any(usable)) usable <- TRUE
    minMaxFn <- if(nodes[usable, NODE_CLASS[1] == nodeClasses$SPLIT]) min else max
    refJxnKey <- nodes[ 
        usable,
        .(minMaxClipLen = minMaxFn(CLIP_LEN)), 
        by = jxnKey
    ][
      order(-minMaxClipLen)
    ][ 
       1,
       jxnKey
    ]
    refMol <- jxnMols[jxnKey == refJxnKey]
    refNodes <- nodes[jxnKey == refJxnKey] # exactly two nodes that help build the call set 
    if(refMol[1, JXN_TYPE] == junctionTypes$UNKNOWN) return(NULL)

    #########################
    return(refMol)


    # add any outer clips that match inner clips at the reference nodes (whether reference in a split or a gap)
    # TODO: is this correct for both svCapture and svWGS? probably need to track IS_COLLATED
    isInversion <- refMol[, side1 == side2] 
    flipNode    <- if(refMol[, side1 == 'R']) 1L else 2L
    for(i in 1:2){
        if(refNodes[i, CLIP_LEN == 0]) next # unclipped gap nodes aren't informative for outer clip matching

        # TODO: need to mount all samples for fast reading of clips, do this for sample in SAMPLES
        outerClipNodes <- getNodes('outer_clips', refNodes[i, NODE], unpackNodeNames = TRUE)
        if(is.null(outerClipNodes)) next
        outerClipNodes[, ':='(
            jxnKey = paste(MOL_ID, JXN_N, sep = ":"),            
            NODE_N = i,
            N_COLLAPSED = 1
        )]

        # TODO: simplify this check - only check duplication at the outer clip end
        # will tend to artifactually collapse some truly different molecules but effect on outcomes is negligible
        # if(nrow(outerClipNodes) > 1) outerClipNodes <- purgeDuplicateMolecules(outerClipNodes)

        # flip inversion outer clips as needed to also put them on the canonical strand
        if(isInversion && flipNode == i){        
            outerClipNodes[, ':='(
                SEQ      = rc(SEQ),
                CLIP_SEQ = rc(CLIP_SEQ),
                CIGAR    = rc_cigar(CIGAR)
            )]
        } 

        # merge inner and outer node evidence
        nodes <- rbind(nodes, outerClipNodes)
    }

#-------------------------------------------------------------------------------------
    # attempt to merge gap+clip-only junctions by using clipped nodes from different molecules
    # could be inner or outer clips of any length, but each partner must be clipped to correctly set 'pos'
#-------------------------------------------------------------------------------------
    if(refNodes[1, NODE_CLASS == nodeClasses$SPLIT]) { 
        call$JXN_SEQ <- refNodes[1, SEQ] 
    } else if(nodes[, length(unique(MOL_ID)) > 1]){ # only relevant when there are multiple molecules but no splits
        
        # get the two nodes with the greatest potential for overlap with the other side
        # NB: merge nodes can include outer clips
        mergeNodes <- do.call('rbind', lapply(1:2, function(i){
            d <- nodes[NODE_N == i]
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
                nodeN <- nodes[i, NODE_N]
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
