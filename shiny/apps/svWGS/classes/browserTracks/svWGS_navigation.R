# svWGS track navigation

svWGS_navTable_display <- function(jxns){
    req(jxns)
    jxns[, .(
        svId = SV_ID,
        samples = trimws(gsub(",", " ", samples)),
        type = svx_jxnTypes[edgeType, name],
        size = SV_SIZE,
        insertSize = -MICROHOM_LEN,
        nSamples,
        nInstances,
        rStart = paste0(cChrom1, ":", cRefPos1, ifelse(cStrand1 == 1, "+", "-")),
        rEnd   = paste0(cChrom2, ":", cRefPos2, ifelse(cStrand2 == 1, "+", "-"))
    )]   
}
