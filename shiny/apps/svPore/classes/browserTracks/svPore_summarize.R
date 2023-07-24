# svPore summary of all filtered junctions
svPore_summarizeJunctions <- function(jxns, track, layout){
    req(jxns)
    startSpinner(session, message = "creating distributions")
    jxns[, alnBaseQual := as.integer(alnBaseQual)]
    jxns[, alnSize_kb  := as.integer(alnSize / 1000)]
    svx_summarizeJunctions(jxns, track, layout, list(
        jxns[, .N, keyby = .(mapQ)],
        jxns[, .N, keyby = .(alnSize_kb)]
    ))
    # gapCompressedIdentity = max(gapCompressedIdentity),
}
