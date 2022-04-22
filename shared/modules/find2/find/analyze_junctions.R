#-------------------------------------------------------------------------------------
# attempt to merge gap+clip-only junctions by using clipped nodes from different molecules
# could be inner or outer clips of any length, but each partner must be clipped to correctly set 'pos'
#-------------------------------------------------------------------------------------
mergeGapJunctions <- function(svIdx){

    # process only relevant when there are multiple independent molecules but no splits
    jxnMols <- jxnMols[svIdx] 
    if(
        jxnMols[, any(NODE_CLASS == nodeClasses$SPLIT)] || 
        jxnMols[, length(unique(jxnKey)) == 1]
    ){
        return(jxnMols)
    }   

    # get the two nodes with the greatest potential for overlap with the other side
    # NB: merge nodes can include outer clips
    refMol <- jxnMols[IS_REFERENCE == 1]
    mergeMols <- rbind(
        jxnMols[ 
            if(refMol[, side1] == 'L') 
                jxnMols[, which.max(ifelse(NODE_1 == "*", NA, pos1 + CLIP_LEN_1))]
            else 
                jxnMols[, which.min(ifelse(NODE_1 == "*", NA, pos1 - CLIP_LEN_1))] 
        ],
        jxnMols[
            if(refMol[, side2] == 'L') 
                jxnMols[, which.max(ifelse(NODE_2 == "*", NA, pos2 + CLIP_LEN_2))]
            else 
                jxnMols[, which.min(ifelse(NODE_2 == "*", NA, pos2 - CLIP_LEN_2))]            
        ]
    )
    if(nrow(mergeMols) != 2) return(jxnMols)
    mergeNodes <- rbind(
        mergeMols[1, .(
            jxnKey   = jxnKey,
            CLIP_LEN = CLIP_LEN_1,
            SEQ      = SEQ_1
        )],
        mergeMols[2, .(
            jxnKey   = jxnKey,
            CLIP_LEN = CLIP_LEN_2,
            SEQ      = SEQ_2
        )]
    )

    # both nodes must be clipped, so that node positions characterize a junction
    # merge nodes must be from different source molecules (we already tried to merge within the same molecule)
    if(mergeNodes[, any(CLIP_LEN == 0)] || 
       mergeNodes[1, jxnKey] == mergeNodes[2, jxnKey]) return(jxnMols)

    # find the best overlap score over all possible merge registers
    maxOverlap <- mergeNodes[, min(nchar(SEQ))]
    if(maxOverlap < env$MIN_MERGE_OVERLAP) return(jxnMols)
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
            MERGE_LEN = overlapLengths[iAtMax]
        )]               
    }                
    jxnMols
}
#-------------------------------------------------------------------------------------
# describe a consensus SV based on sequence data from a set of molecules discovered as evidence
#-------------------------------------------------------------------------------------
characterizeSvJunction <- function(svIdx){
    jxnMols <- jxnMols[svIdx]
    refMolI <- jxnMols[, .I[IS_REFERENCE == 1]]
    refMol  <- jxnMols[refMolI]
    gapSplitIs <- jxnMols[, .I[NODE_CLASS != nodeClasses$OUTER_CLIP]]

    # characterize the SV junction when possible
    jxnMols[, ':='(
        MICROHOM_LEN = 0L,
        JXN_BASES = "*"
    )]
    if(refMol[, JXN_SEQ != "*"]){
        
        # positions within the complete junction sequence as bounded by refNode1 and refNode2 SEQ limits
        seq1Len  <- refMol[, nchar(SEQ_1)]
        readPos1 <- seq1Len - refMol[, CLIP_LEN_1]
        readPos2 <- refMol[, 
            if(MERGE_LEN == 0) CLIP_LEN_2             # a sequenced junction
            else seq1Len - MERGE_LEN + CLIP_LEN_2 + 1 # a reconstructed junction
        ]

        # parse microhomology vs. inserted bases
        jxnMols[refMolI, MICROHOM_LEN := readPos1 - readPos2 + 1]  
        jxnMols[refMolI, 
            JXN_BASES := if(MICROHOM_LEN > 0) substr(JXN_SEQ, readPos2, readPos1)         # microhomology
                    else if(MICROHOM_LEN < 0) substr(JXN_SEQ, readPos1 + 1, readPos2 - 1) # inserted bases
                    else "" # a blunt joint
        ]   
    }

    # aggregate SV-level information
    jxnMols[,.(
        MAPQ_1          = paste(MAPQ_1, collapse = ","),
        UMI_1           = paste(UMI_1, collapse = ","),
        #-------------
        N_TOTAL         = .N,
        N_GAPS          = sum(NODE_CLASS == nodeClasses$GAP),
        N_SPLITS        = sum(NODE_CLASS == nodeClasses$SPLIT),
        N_OUTER_CLIPS   = sum(NODE_CLASS == nodeClasses$CLIP),
        #-------------
        JXN_TYPE        = JXN_TYPE[refMolI],
        #-------------
        N_DUPLEX        = sum(IS_DUPLEX), # including clips
        N_DUPLEX_GS     = sum(IS_DUPLEX[gapSplitIs]), # splits + gaps only
        STRAND_COUNT    = sum(STRAND_COUNT1 + STRAND_COUNT2),
        STRAND_COUNT_GS = sum(STRAND_COUNT1[gapSplitIs] + STRAND_COUNT2[gapSplitIs]),
        STRAND_COUNT1   = sum(STRAND_COUNT1),
        STRAND_COUNT2   = sum(STRAND_COUNT2),
        TARGET_CLASS    = TARGET_CLASS[refMolI],    
        SHARED_PROPER   = mean(SHARED_PROPER),
        SHARED_PROPER_GS= mean(SHARED_PROPER[gapSplitIs]),
        #------------- 
        SAMPLES         = paste(unique(SAMPLE), collapse = ","),
        N_SAMPLES       = length(unique(SAMPLE)),     
        #------------- 
        MAPQ_2          = paste(MAPQ_2, collapse = ","),
        UMI_2           = paste(UMI_2, collapse = ","),
        #------------- 
        CHROM_1         = revChromIndex[[chrom1[refMolI]]],
        SIDE_1          = side1[refMolI],
        POS_1           = pos1[refMolI],
        CHROM_2         = revChromIndex[[chrom2[refMolI]]],
        SIDE_2          = side2[refMolI],
        POS_2           = pos2[refMolI],
        #------------- 
        JUNCTION_NAME   = jxnName[refMolI],
        JUNCTION_NAMES  = paste(unique(jxnName), collapse = "::"), # includes JUNCTION_NAME
        #------------- 
        SV_ID           = svIndex[refMolI],
        #------------- 
        N_AMBIGUOUS     = sum(AMBIGUOUS),
        N_DOWNSAMPLED   = sum(DOWNSAMPLED),
        N_COLLAPSED     = sum(N_COLLAPSED > 1),
        #------------- 
        JXN_SEQ         = JXN_SEQ[refMolI],
        MERGE_LEN       = MERGE_LEN[refMolI],
        #-------------
        MICROHOM_LEN    = MICROHOM_LEN[refMolI],
        JXN_BASES       = JXN_BASES[refMolI],
        #-------------
        SV_SIZE         = if(JXN_TYPE[refMolI] == "T") 0 else abs(pos2[refMolI] - pos1[refMolI]),
        #-------------
        GEN_REF_1       = getRefSeq_padded(revChromIndex[[chrom1[refMolI]]], pos1[refMolI], faidx_padding),
        GEN_REF_2       = getRefSeq_padded(revChromIndex[[chrom2[refMolI]]], pos2[refMolI], faidx_padding),
        #-------------
        TARGET_REGION = paste(
            sort(unique(c(
                getTargetRegionName(chrom1[refMolI], pos1[refMolI]),
                getTargetRegionName(chrom2[refMolI], pos2[refMolI])
            ))),
            collapse = ","
        ),
        TARGET_POS_1    = getTargetRegionI(chrom1[refMolI], pos1[refMolI]),
        TARGET_POS_2    = getTargetRegionI(chrom2[refMolI], pos2[refMolI]) 
    )]
}
#=====================================================================================
