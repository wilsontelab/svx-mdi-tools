#----------------------------------------------------------------------
# handle unique junction loading and filtering
#----------------------------------------------------------------------

# app-specific junction loading
svPore_loadJunctions <- function(sourceId){
    jxns <- readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
    jxns[, ":="(
        nSequenced = nInstances, # for nanopore, all junctions were sequenced (there are not junction in pair gaps)
        flankLength = alnSize,
        nLinkedJunctions = segmentClusterN
    )]
    jxns
#  [1] "clusterN"                "edgeType"
#  [3] "eventSize"               "node1"
#  [5] "node2"                   "cChromIndex1"
#  [7] "cChromIndex2"            "cStrand1"
#  [9] "cStrand2"                "cRefPos1"
# [11] "cRefPos2"                "insertSize"
# [13] "mapQ"                    "gapCompressedIdentity"
# [15] "baseQual"                "alnBaseQual"
# [17] "alnSize"                 "samples"
# [19] "nSamples"                "nInstances"
# [21] "nCanonical"              "nNonCanonical"
# [23] "HCT_APH_untargeted"      "HCT_APH_untargeted_redo"
# [25] "segmentClusterN"
}

# get a single junction to expand
svPore_getJunction <- function(x){
    svx_loadJunctions(x$sourceId, svPore_loadJunctions)[clusterN == x$clusterN]
}
