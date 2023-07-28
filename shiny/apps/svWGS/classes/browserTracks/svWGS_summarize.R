# svWGS summary of all filtered junctions
svWGS_summarizeJunctions <- function(jxns, track, layout){
    req(jxns)
    startSpinner(session, message = "creating distributions")
    jxns[, insertSize := -MICROHOM_LEN]
    jxns[, nSequencedJxn := N_SPLITS]
    jxns[, mergeLen := MERGE_LEN]
    svx_summarizeJunctions(jxns, track, layout, list(
        jxns[, .N, keyby = .(mapQ)],
        jxns[, .N, keyby = .(nSequencedJxn)],
        jxns[, .N, keyby = .(mergeLen)]
    ))
}
