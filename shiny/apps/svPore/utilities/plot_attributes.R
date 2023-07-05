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
    longName = c(
        "Deletion",
        "Insertion",
        "Duplication",
        "Inversion",
        "Translocation",
        "?"
    ), 
    legend = c(
        TRUE,
        FALSE,
        TRUE,
        TRUE,
        TRUE,
        FALSE
    ),
    order = 1:6,
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
svPore$setJunctionPointColors <- function(jc, track){
    jc[, color := switch(
        getBrowserTrackSetting(track, "Points", "Color_By", "edgeType"),
        edgeType = svPore$jxnTypes[edgeType, color],
        isSingleton    = ifelse(nInstances == 1, CONSTANTS$plotlyColors$red, CONSTANTS$plotlyColors$black),
        isSingleSample = ifelse(nSamples == 1,   CONSTANTS$plotlyColors$red, CONSTANTS$plotlyColors$black),
        "black" # will be corrected later when Color_By == "sample"
    )]
    jc
}
svPore$setJunctionPointSizes <- function(jc, track){
    cex1      <- getBrowserTrackSetting(track, "Points", "Point_Size", 0.25)
    cexFactor <- getBrowserTrackSetting(track, "Points", "Point_Scale_Factor", 0.25)
    jc[, cex := switch(
        getBrowserTrackSetting(track,"Points","Size_By","fixed"),
        fixed = cex1,
        nSamples   = cex1 * nSamples   * cexFactor,
        nInstances = cex1 * nInstances * cexFactor,
        nSampleInstances = cex1 * nSampleInstances * cexFactor,
        0.25
    )]
    jc
}
svPore$getSampleColors <- function(sourcesToPlot){
    allSamples <- unique(unlist(lapply(sourcesToPlot, function(x) x$Sample_ID)))
    sampleCols <- as.list(1:length(allSamples))
    names(sampleCols) <- paste0(",", allSamples, ",") 
    sampleCols   
}
svPore$colorBySample <- function(jc, sourcesToPlot){
    sampleCols <- svPore$getSampleColors(sourcesToPlot)
    jc[, color := ifelse(
        nSamples > 1, 
        "black", 
        unlist(CONSTANTS$plotlyColors[unlist(sampleCols[samples])])
    )]
}

svPore$junctionClusterLegend <- function(track, coord, ylim, sourcesToPlot){
    x <- switch(
        getBrowserTrackSetting(track, "Points", "Color_By", "edgeType"),
        edgeType = {
            jt <- svPore$jxnTypes[order(order)]
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
            sampleCols <- svPore$getSampleColors(sourcesToPlot)            
            list(
                legend = c(gsub(",", "", names(sampleCols)), ">1 sample"),
                color  = c(unlist(CONSTANTS$plotlyColors[unlist(sampleCols)]), CONSTANTS$plotlyColors$black)
            )
        }
    )
    trackLegend(
        track, coord, ylim, 
        legend = x$legend, pch = 19, col = x$color
    )
}
