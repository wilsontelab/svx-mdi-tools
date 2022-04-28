#-------------------------------------------------------------------------------------
# examine a set of SV-nominating molecules in a specific sample for nearly-identical molecule spans
# aggregate presumed duplicate molecules to prevent them from appearing as falsely independent evidence
#-------------------------------------------------------------------------------------
purgeDuplicateMolecules <- function(smpSvIdx){
    purgeDuplicateMolecules_(jxnMols[smpSvIdx])
}
purgeDuplicateMolecules_ <- function(jxnMols){

    # check if singleton-molecule call
    if(nrow(jxnMols) == 1) return(jxnMols)

    # randomly downsample junctions with extreme coverage for processing efficiency
    # this will still allow robust SV calling, but does impair accurate coverage determination
    if(nrow(jxnMols) > env$PURGE_LIMIT){
        jxnMols <- jxnMols[sample(.N, env$PURGE_LIMIT)]
        jxnMols[, DOWNSAMPLED := 1]
    }

    # use the euclidean distance between the endpoints of each pair of molecules to find collisions
    setkey(jxnMols, jxnKey)
    pairs <- as.data.table(t(combn(jxnMols[, jxnKey], 2)), stringsAsFactors = FALSE)
    setnames(pairs, c('jxnKey1', 'jxnKey2'))        
    pairs <- cbind(
        pairs, 
        jxnMols[pairs[, jxnKey1], .(OUT_POS1_1 = OUT_POS1, OUT_POS2_1 = OUT_POS2)], 
        jxnMols[pairs[, jxnKey2], .(OUT_POS1_2 = OUT_POS1, OUT_POS2_2 = OUT_POS2)]
    )
    collisionIs <- pairs[, .I[sqrt((OUT_POS1_1 - OUT_POS1_2) ** 2 + 
                                   (OUT_POS2_1 - OUT_POS2_2) ** 2) <= env$PURGE_DISTANCE]]

    # check if no collisions
    if(length(collisionIs) == 0) return(jxnMols)

    # for all colliding pairs, keep the molecule with the best read evidence and add other molecules to its counts
    remappedJxnKeys <- list()
    rejectedMolIds <- integer()
    for(i in collisionIs){ 

        # determine which molecule of the colliding pair to keep
        pair <- jxnMols[pairs[i, c(jxnKey1, jxnKey2)]][order(-NODE_CLASS, -(STRAND_COUNT1 + STRAND_COUNT2), -(MAPQ_1 + MAPQ_2))]
        bestJxnKey  <- pair[1, jxnKey]
        worstJxnKey <- pair[2, jxnKey]

        # ascend up any prior chain to the best of a set of colliding molecules
        while(!is.null(remappedJxnKeys[[bestJxnKey]])) bestJxnKey <- remappedJxnKeys[[bestJxnKey]]
        remappedJxnKeys[[worstJxnKey]] <- bestJxnKey # the molecule that holds a rejected molecule's counts

        # add the worst strand counts to the best molecule and count how often we did this
        worstCount1 <- pair[2, STRAND_COUNT1]
        worstCount2 <- pair[2, STRAND_COUNT2]
        jxnMols[bestJxnKey, ':='(
            STRAND_COUNT1 = STRAND_COUNT1 + worstCount1,
            STRAND_COUNT2 = STRAND_COUNT2 + worstCount2,
            N_COLLAPSED   = N_COLLAPSED + 1L
        )]

        # add to the list of rejected molecules
        rejectedMolIds <- append(rejectedMolIds, pair[2, MOL_ID])
    }   

    # finally, completely purge the presumed duplicate _molecules_ from the nodes list
    jxnMols <- jxnMols[!(MOL_ID %in% rejectedMolIds)]
    jxnMols
}
