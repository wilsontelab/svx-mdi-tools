#-------------------------------------------------------------------------------------
# select one molecule as a reference for characterizing the junction
#   prefer splits, which can always completely sequence a junction  
#   for splits, prefer longest clip on the shorter-clip side (promotes long molecules with central junctions)
#   for gaps, prefer molecules with the longest clips on either side
#-------------------------------------------------------------------------------------
assignReferenceMolecule <- function(svIdx){
    jxnMols <- jxnMols[svIndex == svIdx] 
    if(nrow(jxnMols) == 1){
        jxnMols[, IS_REFERENCE := 1]
        return(jxnMols)
    }   
    usable <- jxnMols[, NODE_CLASS == nodeClasses$SPLIT]
    if(!any(usable)) usable <- TRUE
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
# print clipped nodes from reference molecules to support outer clip recovery
printReferenceNodes <- function(){
    refNodes <- do.call(rbind, lapply(1:2, function(i){
        node    <- paste('NODE',     i, sep = "_")
        clipLen <- paste('CLIP_LEN', i, sep = "_")
        x <- jxnMols[refMols, .SD, .SDcols = c('svIndex', node, clipLen)]
        setnames(x, c('svIndex', 'NODE', 'CLIP_LEN'))
        x[, NODE_N := i]
        x
    }))
    refNodesFile <- paste(env$FIND_PREFIX, 'clipped_reference_nodes', 'txt', sep = ".")
    write.table(
        refNodes[CLIP_LEN > 0, .(svIndex, NODE, NODE_N)],  # don't look for clip matches to unclipped gap inner nodes
        file = refNodesFile, 
        quote = FALSE, 
        sep = "\t",    
        row.names = FALSE,   
        col.names = FALSE
    )
    list(
        file = refNodesFile,
        table = refNodes # all ref nodes, whether clipped or not
    )
}
