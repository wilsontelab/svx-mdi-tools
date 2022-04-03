#-------------------------------------------------------------------------------------
# select one molecule as a reference for characterizing the junction
#   prefer splits, which can always completely sequence a junction  
#   for splits, prefer longest clip on the shorter-clip side (promotes long molecules with central junctions)
#   for gaps, prefer molecules with the longest clips on either side
#-------------------------------------------------------------------------------------
assignReferenceMolecule <- function(svIdx){
    jxnMols <- jxnMols[svIdx] 
    if(nrow(jxnMols) == 1){
        jxnMols[, IS_REFERENCE := 1]
        return(jxnMols)
    }   
    usable <- jxnMols[, NODE_CLASS == nodeClasses$SPLIT]
    if(!any(usable)) usable <- TRUE # i.e., all gap molecules continue if no splits
    nodes <- rbind(
        jxnMols[usable, .(jxnKey = jxnKey, CLIP_LEN = CLIP_LEN_1)],
        jxnMols[usable, .(jxnKey = jxnKey, CLIP_LEN = CLIP_LEN_2)]
    )
    minMaxFn <- if(jxnMols[usable, NODE_CLASS[1] == nodeClasses$SPLIT]) min else max
    refJxnKey <- nodes[,
        .(minMaxClipLen = minMaxFn(CLIP_LEN)), 
        by = jxnKey
    ][
      order(-minMaxClipLen)
    ][ 
       1,
       jxnKey
    ]
    jxnMols[, IS_REFERENCE := ifelse(jxnKey == refJxnKey, 1, 0)]
    jxnMols
}

# print clipped nodes from reference molecules for subsequent outer clip recovery
printReferenceNodes <- function(){
    refNodes <- do.call(rbind, lapply(1:2, function(i){
        node      <- paste('NODE',     i, sep = "_")
        clipLen   <- paste('CLIP_LEN', i, sep = "_")
        otherNode <- paste('NODE', i %% 2 + 1, sep = "_")
        x <- jxnMols[IS_REFERENCE == 1, .SD, .SDcols = c('svIndex', node, clipLen, otherNode)]
        setnames(x, c('svIndex', 'NODE', 'CLIP_LEN', 'OTHER_NODE'))
        x[, NODE_N := i]
        x
    }))
    refNodesFile <- paste(env$FIND_PREFIX, 'clipped_reference_nodes', 'txt', sep = ".")
    write.table( # don't look for outer clip matches to unclipped gap inner nodes
        refNodes[CLIP_LEN > 0, .(svIndex, NODE, NODE_N, OTHER_NODE)],
        file = refNodesFile, 
        quote = FALSE, 
        sep = "\t",    
        row.names = FALSE,   
        col.names = FALSE
    )
    refNodesFile
}
