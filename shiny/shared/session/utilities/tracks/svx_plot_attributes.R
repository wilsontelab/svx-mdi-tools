# junction plot attributes
svx_setJunctionPointColors <- function(jxns, track, family = "Points"){
    jxns[, color := switch(
        getBrowserTrackSetting(track, family, "Color_By", "edgeType"),
        edgeType = svx_jxnType_codeToX(edgeType, "color"),
        isSingleton    = ifelse(nInstances == 1, CONSTANTS$plotlyColors$red, CONSTANTS$plotlyColors$black),
        isSingleSample = ifelse(nSamples == 1,   CONSTANTS$plotlyColors$red, CONSTANTS$plotlyColors$black),
        "black" # will be corrected later when Color_By == "sample"
    )]
    jxns
}
svx_setJunctionPointSizes <- function(jxns, track, family = "Points"){
    cex1      <- getBrowserTrackSetting(track, family, "Point_Size", 0.5)
    cexFactor <- getBrowserTrackSetting(track, family, "Point_Scale_Factor", 0.25)
    jxns[, cex := switch(
        getBrowserTrackSetting(track,family,"Size_By","fixed"),
        fixed = cex1,
        nSamples   = cex1 * nSamples   * cexFactor,
        nInstances = cex1 * nInstances * cexFactor,
        nSampleInstances = cex1 * nSampleInstances * cexFactor,
        0.25
    )]
    jxns
}

# junctions legend
svx_junctionsLegend <- function(track, coord, ylim, selectedTargets, 
                                isMultiSample = TRUE, jxns = NULL, family = "Points", sampleNameFn = NULL){
    x <- switch(
        getBrowserTrackSetting(track, family, "Color_By", "edgeType"),
        edgeType = {
            jt <- svx_jxnTypes[order(order)]
            list(
                legend = jt$longName[jt$legend],
                color  = jt$color[jt$legend]
            )
        },
        isSingleton = list(
            legend = c(">1 molecule","1 molecule"),
            color  = c(CONSTANTS$plotlyColors$black, CONSTANTS$plotlyColors$red)
        ),
        isSingleSample = list(
            legend = c(">1 sample","1 sample"),
            color  = c(CONSTANTS$plotlyColors$black, CONSTANTS$plotlyColors$red)
        ),
        sample = {
            sampleCols <- getColorsBySelectedSample(selectedTargets, isMultiSample, jxns) 
            sampleNames <- names(sampleCols) 
            if(!is.null(sampleNameFn)) sampleNames <- sampleNameFn(sampleNames, selectedTargets, jxns)
            list(
                legend = c(gsub(",", "", sampleNames), ">1 sample"),
                color  = c(unlist(CONSTANTS$plotlyColors[unlist(sampleCols)]), CONSTANTS$plotlyColors$black)
            )
        }
    )
    trackLegend(
        track, coord, ylim, 
        legend = x$legend, pch = 19, pt.cex = 1.25, col = x$color
    )
}
