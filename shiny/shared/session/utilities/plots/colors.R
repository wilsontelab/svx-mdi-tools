#----------------------------------------------------------------------
# set SV plot point colors
#----------------------------------------------------------------------
getSvPointColors <- function(filteredSvs, settings, sampleSelector, isCapture = FALSE){
    reactive({
        svs <- filteredSvs()
        req(svs)
        ps <- settings$Plot_Settings()
        blue <- list(
            colors = CONSTANTS$plotlyColors$blue, 
            color = NULL,
            label = NULL
        )
        switch(ps$Point_Color$value,
            type = {
                svFilters <- settings$SV_Filters()
                list(
                    colors = svs[, svx_jxnType_codeToX(edgeType, "color")],
                    color = svx_jxnTypes[legend == TRUE & name %in% svFilters$SV_Type$value, color],
                    label = svx_jxnTypes[legend == TRUE & name %in% svFilters$SV_Type$value, name]
                )                
            },
            sample = {
                assignments <- sampleSelector$selectedAssignments()
                svs_1 <- svs[N_SAMPLES == 1]
                svs_1[, uniqueId := unlist(PROJECT_SAMPLES)]
                svs_1 <- merge(
                    svs_1, 
                    assignments[, .(uniqueId, Category1, Category2)],
                    all.x = TRUE,
                    by = "uniqueId"
                )
                samples <- svs_1[, sort(unique(SAMPLES)), keyby = c("Category1", "Category2")][[3]]
                sampleColors <- unlist(CONSTANTS$plotlyColors[seq_along(samples)])
                names(sampleColors) <- samples     
                sampleLabels <- sapply(samples, function(sample){
                    id <- assignments[Sample_ID == sample, uniqueId]
                    getSampleNames(sampleUniqueIds = id)
                })
                list(
                    colors = svs[, ifelse(
                        N_SAMPLES == 1, 
                        sampleColors[SAMPLES], 
                        CONSTANTS$plotlyColors$grey
                    )],
                    color = sampleColors,
                    label = sampleLabels
                )
            },
            target = {
                if(isCapture){
                    targets <- svs[!grepl(',', TARGET_REGION), sort(unique(TARGET_REGION))]
                    targetColors <- unlist(CONSTANTS$plotlyColors[seq_along(targets)])
                    names(targetColors) <- targets
                    list(
                        colors = svs[, ifelse(
                            grepl(',', TARGET_REGION), 
                            CONSTANTS$plotlyColors$grey, 
                            targetColors[TARGET_REGION])
                        ],
                        color = targetColors,
                        label = targets
                    )
                } else blue
            }, 
            duplex = list(
                colors = svs[, ifelse(
                    N_DUPLEX_GS > 0, 
                    CONSTANTS$plotlyColors$green,  
                    CONSTANTS$plotlyColors$red
                )],
                color = c(
                    CONSTANTS$plotlyColors$green,  
                    CONSTANTS$plotlyColors$red
                    ),
                label = c(
                    "Duplex", 
                    "Single Strand"
                )
            ),
            blue = blue
        )
    })
}
pointColorLegend <- function(stepSettings, plotSettings, svPointColors, exclude = character()){
    if(is.null(svPointColors$label)) return()
    is <- which(!(svPointColors$label %in% exclude))
    legend(
        plotSettings$get('Plot_Frame', 'Legend_Placement'),
        svPointColors$label[is],
        pch = 20, 
        pt.cex = stepSettings$Point_Size$value * 1.25,
        col = svPointColors$color[is]
    )
}
