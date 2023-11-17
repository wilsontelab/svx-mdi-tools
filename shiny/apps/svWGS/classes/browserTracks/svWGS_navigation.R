# svWGS track navigation

svWGS_navTable_display <- function(jxns){
    req(jxns)
    jxns[, .(
        svId = SV_ID,
        samples = trimws(gsub(",", " ", samples)),
        type = svx_jxnType_codeToX(edgeType, "name"),
        cnc = round(maxCNC, 1),
        size = SV_SIZE,
        insertSize,
        nSamples,
        nInstances,
        nSequenced,
        flankLength,
        nLinkedJxns = nLinkedJunctions,
        rStart = paste0(cChrom1, ":", cRefPos1, ifelse(cStrand1 == 1, "+", "-")),
        rEnd   = paste0(cChrom2, ":", cRefPos2, ifelse(cStrand2 == 1, "+", "-"))
    )]   
}
