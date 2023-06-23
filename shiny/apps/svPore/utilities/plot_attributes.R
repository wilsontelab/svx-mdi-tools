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
svPore$setJunctionPointColors <- function(j, settings){
    j[, color := switch(
        settings$get("Points", "Color_By"),
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
svPore$setJunctionPointSizes <- function(j, settings){
    cex1 <- settings$get("Points", "Point_Size")
    j[, cex := switch(
        settings$get("Points", "Size_By"),
        fixed = cex1,
        nSamples   = cex1 * nSamples,
        nInstances = cex1 * nInstances,
        nSampleInstances = cex1 * nSampleInstances,
        0.25
    )]
    j
}