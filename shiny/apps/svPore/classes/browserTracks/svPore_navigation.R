# svPore track navigation
svx_navTable_display_app <- function(jxns){
    req(jxns)
    jxns[, .(
        cluster = clusterN,
        samples = trimws(gsub(",", " ", samples)),
        type = svx_jxnTypes[edgeType, name],
        size = eventSize,
        insertSize,
        nSamples,
        nInstances,
        rStart = paste0(cChrom1, ":", cRefPos1, ifelse(cStrand1 == 1, "+", "-")),
        rEnd   = paste0(cChrom2, ":", cRefPos2, ifelse(cStrand2 == 1, "+", "-"))
    )]   
}
