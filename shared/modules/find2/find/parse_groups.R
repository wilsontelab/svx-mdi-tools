#-------------------------------------------------------------------------------------
# process continuity groups to find SV matches between junction molecules
# these steps are executed without respect to source sample
#-------------------------------------------------------------------------------------

# label and return the junction molecules that are evidence of a single SV call
commitJxnMols <- function(jxnMols, svJ){

    # apply the SV coverage filter
    if(nrow(jxnMols) < env$MIN_COVERAGE) return(NULL)

    # enforce the SV-level MAPQ filter (generally higher than the original molecule-level filter)
    maxMapQ1 <- jxnMols[, max(MAPQ_1)]
    maxMapQ2 <- jxnMols[, max(MAPQ_2)]
    if(jxnMols[, # filter is applied to the _best_ gap or split in the molecule set
        !(maxMapQ1 >= env$MIN_MAPQ_BOTH && maxMapQ2 >= env$MIN_MAPQ_BOTH) ||
        !(maxMapQ1 >= env$MIN_MAPQ_ONE  || maxMapQ2 >= env$MIN_MAPQ_ONE)    
    ]) return(NULL)

    # all filters passed, set SV ids and return
    jxnMols[, ':='(
        svIndex = paste(groupIndex, svJ, sep = ":"),
        sampleSvIndex = paste(SAMPLES[SAMPLE], groupIndex, svJ, sep = ":")
    )] 
    jxnMols
}

# compare a gap junction in one molecule to a sequenced junction in another molecule
checkHalfSequenced <- function(delta, side) { # delta must be calculated as SPLIT - GAP positions
    if(side == "L") delta >= -env$PURGE_DISTANCE && delta <= MAX_MAX_TLEN
               else delta >= -MAX_MAX_TLEN && delta <= env$PURGE_DISTANCE  
}
# determine if edges in different molecules appear to be the same SV on one side of the junction
isNodeMatch <- function(IS_SPLIT_J1, IS_SPLIT_J2, POS_J1, POS_J2, side){
         if( IS_SPLIT_J1 &&  IS_SPLIT_J2) abs(POS_J1 - POS_J2) <= env$PURGE_DISTANCE
    else if(!IS_SPLIT_J1 && !IS_SPLIT_J2) abs(POS_J1 - POS_J2) <= MAX_MAX_TLEN # TODO: use sample-specific MAX_TLEN?
    else if(IS_SPLIT_J1) checkHalfSequenced(POS_J1 - POS_J2, side)
                    else checkHalfSequenced(POS_J2 - POS_J1, side)
}
# determine which jxnMols match a seed junction
addMatches <- function(matches, jxnTypes, seedJxnType, side1, side2){
    SEED_IS_SPLIT <- jxnTypes[seedJxnType, NODE_CLASS == nodeClasses$SPLIT]
    IS_SPLIT      <- jxnTypes[           , NODE_CLASS == nodeClasses$SPLIT] 
    SEED_POS1     <- jxnTypes[seedJxnType, pos1] # node inner positions, not alignment positions
         POS1     <- jxnTypes[           , pos1]
    SEED_POS2     <- jxnTypes[seedJxnType, pos2]
         POS2     <- jxnTypes[           , pos2]
    cbind(matches, 
        (SEED_IS_SPLIT | !IS_SPLIT) & # never match splits to gap seeds
        mapply( # position collison on side1
            isNodeMatch, 
            SEED_IS_SPLIT, 
            IS_SPLIT, 
            SEED_POS1, 
            POS1,
            side1
        ) & 
        mapply( # position collison on side2
            isNodeMatch, 
            SEED_IS_SPLIT, 
            IS_SPLIT, 
            SEED_POS2, 
            POS2,
            side2
        ) & 
        abs(SEED_POS1 - POS1) + abs(SEED_POS2 - POS2) <= MAX_MAX_TLEN # ensure inferred gap molecule size is not be too large
    )
}

# determine which junction molecules in a continuity group appear to be the same SV
parseContinuityGroup <- function(grpIdx){

    # check if singleton-molecule call or a single-junction set, e.g., multiple identical splits
    jxnMols <- jxnMols[groupIndex == grpIdx]
    if(nrow(jxnMols) == 1 ||
       jxnMols[, length(unique(jxnName)) == 1]) return(commitJxnMols(jxnMols, 1))

    # collapse junction molecules across multiple samples to junction types = NODE_CLASS:NODE_1:NODE_2
    jxnTypes <- jxnMols[, .(
        jxnType      = paste(NODE_CLASS[1], jxnName[1], sep=":"),
        jxnKeys      = list(jxnKey), # keep track of all sources molecules with this jxnType        
        pos1         = pos1[1], # all junctions with the same junction name have the same node positions
        pos2         = pos2[1],
        STRAND_COUNT = sum(STRAND_COUNT1 + STRAND_COUNT2), # total _read-pair_ count over all samples
        MOL_COUNT    = .N, # total _molecule_ count over all samples (for uncollated inputs)
        MAPQ         = max(MAPQ_1 + MAPQ_2) # additional tie-breaking based on MAPQ for picking seed junctions
    ), by = c('NODE_CLASS', 'jxnName')]
    setkey(jxnTypes, jxnType)    

    # pick a seed junction, preferring splits > gaps, higher coverage, higher map quality
    seedJxnType <- jxnTypes[order(-NODE_CLASS, -STRAND_COUNT, -MOL_COUNT, -MAPQ)][1, jxnType]

    # mark matches to seed
    matches <- data.table()
    matches <- addMatches(matches, jxnTypes, seedJxnType, # sides are the same throughout a continuity group
                          jxnMols[1, side1], jxnMols[1, side2])

    # check if jxnMols are all one coherent SV group through any combination of splits and gaps
    if(all(matches)) return(commitJxnMols(jxnMols, 1))
    
    # pick the next best seed from the unmatched junction types, iterate until all types match a seed
    matchedJxnTypes <- apply(matches, 1, any)
    while(!all(matchedJxnTypes)){
        seedJxnType <- jxnTypes[matchedJxnTypes == FALSE][order(-NODE_CLASS, -STRAND_COUNT, -MOL_COUNT, -MAPQ)][1, jxnType]
        matches <- addMatches(matches, jxnTypes, seedJxnType, 
                              jxnMols[1, side1], jxnMols[1, side2])
        matchedJxnTypes <- apply(matches, 1, any)
    }

    # analyze the SV groups for ambiguous junction molecules consistent with multiple seeds
    # these arise when:
    #   two junctions are closely spaced, such that a gap molecule could have come from either allele source
    #   an extremely large source molecule failed MAX_TLEN matching, such that it became a false seed
    #   a split node has an alignment position error more extreme than PURGE_DISTANCE, another false seed
    #       this is common when alignment errors at simple repeats generate chains of alignment positions
    
    # first, concatenate split alignment chains into one inferred true junction (splits can never be ambiguous)
    ambiguousJxnTypes <- apply(matches, 1, sum) > 1
    while(jxnTypes[ambiguousJxnTypes, any(NODE_CLASS == nodeClasses$SPLIT)]){  
        row <- which(ambiguousJxnTypes & jxnTypes[, NODE_CLASS == nodeClasses$SPLIT])[1]
        cols <- which(unlist(matches[row]))[1:2] # identify two columns (seed junctions) with a split collision ...
        matches[[cols[1]]] <- matches[[cols[1]]] | matches[[cols[2]]]
        matches[[cols[2]]] <- NULL # ... and merge them
        ambiguousJxnTypes <- apply(matches, 1, sum) > 1
    }

    # commit every (merged) seed junction group as an SV; ambiguous molecules will be committed multiple times
    do.call(rbind, lapply(1:ncol(matches), function(j){
        jxnKeys <- unlist( jxnTypes[matches[[j]], jxnKeys] )
        commitJxnMols(jxnMols[jxnKey %in% jxnKeys], j)
    }))
}
