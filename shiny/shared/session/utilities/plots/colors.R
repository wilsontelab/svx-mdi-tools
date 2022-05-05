#----------------------------------------------------------------------
# set SV plot point colors
#----------------------------------------------------------------------
getSvPointColors <- function(filteredSvs, settings, isCapture = FALSE){
    reactive({
        svs <- filteredSvs()
        req(svs)
        ps <- settings$Plot_Settings()
        setkey(SVX$jxnTypes, code)
        blue <- list(
            colors = CONSTANTS$plotlyColors$blue, 
            color = NULL,
            label = NULL
        )
        switch(ps$Point_Color$value,
            type = {
                svFilters <- settings$SV_Filters()
                list(
                    colors = svs[, SVX$jxnTypes[JXN_TYPE, color]],
                    color = SVX$jxnTypes[!is.na(color) & name %in% svFilters$SV_Type$value, color],
                    label = SVX$jxnTypes[!is.na(color) & name %in% svFilters$SV_Type$value, name]
                )                
            },
            sample = {
                samples <- svs[N_SAMPLES == 1, sort(unique(SAMPLES))]
                sampleColors <- seq_along(samples)
                names(sampleColors) <- samples
                list(
                    colors = svs[, ifelse(
                        N_SAMPLES == 1, 
                        sampleColors[SAMPLES], 
                        CONSTANTS$plotlyColors$grey
                    )],
                    color = sampleColors,
                    label = samples
                )
            },
            target = {
                if(isCapture){
                    targets <- svs[!grepl(',', TARGET_REGION), sort(unique(TARGET_REGION))]
                    targetColors <- seq_along(targets)
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
