#----------------------------------------------------------------------
# types
#----------------------------------------------------------------------
svPore <- list()
svPore$jxnTypes <- data.table(
    code = c(
        "D",
        "I",
        "U",
        "V",
        "T",
        "?"
    ),
    name = c(
        "Del",
        "Ins",
        "Dup",
        "Inv",
        "Trans",
        "?"
    ),
    color = c(
        CONSTANTS$plotlyColors$blue,
        CONSTANTS$plotlyColors$black,
        CONSTANTS$plotlyColors$green,
        CONSTANTS$plotlyColors$red,
        CONSTANTS$plotlyColors$orange,
        NA
    )
)
setkey(svPore$jxnTypes, code)
svPore$setJunctionPointColors <- function(j, track){
    j[, color := switch(
        getBrowserTrackSetting(track, "Points", "Color_By", "edgeType"),
        edgeType = svPore$jxnTypes[edgeType, color],
        isSingleton    = ifelse(nInstances == 1, CONSTANTS$plotlyColors$red, CONSTANTS$plotlyColors$black),
        isSingleSample = ifelse(nSamples == 1,   CONSTANTS$plotlyColors$red, CONSTANTS$plotlyColors$black),
        "black"
        # ,
        # sample = {
        #     factor(sample)

        # }
    )]
    j
}
svPore$setJunctionPointSizes <- function(j, track){
    cex1      <- getBrowserTrackSetting(track, "Points", "Point_Size", 0.25)
    cexFactor <- getBrowserTrackSetting(track, "Points", "Point_Scale_Factor", 0.25)
    j[, cex := switch(
        getBrowserTrackSetting(track,"Points","Size_By","fixed"),
        fixed = cex1,
        nSamples   = cex1 * nSamples   * cexFactor,
        nInstances = cex1 * nInstances * cexFactor,
        nSampleInstances = cex1 * nSampleInstances * cexFactor,
        0.25
    )]
    j
}