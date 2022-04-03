#----------------------------------------------------------
# set plot positions and parameters for SV summary plots
#----------------------------------------------------------
setMoleculeOuterPositions <- function(jxnMols){

    # jxnMols[NODE_CLASS == nodeClasses$OUTER_CLIP, ':='(
    #     OUT_POS1 = ifelse()   getTargetRegionsI(chrom1, OUT_POS1),
    #     OUT_POS2 = getTargetRegionsI(chrom2, OUT_POS2) 
    # )]

    jxnMols[, ':='(
        TARGET_POS_1 = getTargetRegionsI(chrom1, OUT_POS1),
        TARGET_POS_2 = getTargetRegionsI(chrom2, OUT_POS2) 
    )]
}
