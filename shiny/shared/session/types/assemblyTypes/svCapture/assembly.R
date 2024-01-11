# ----------------------------------------------------------------------
# top-level support for svCapture assembly loading and parsing
# ----------------------------------------------------------------------
svCapture_getDenominatorColumn <- function(settings){
    x <- "kbp_on_target"
    if(isTruthy(settings$get("SV_Capture","Enforce_Min_Family_Size"))) x <- paste(x, "filtered", sep = "_")
    if(isTruthy(settings$get("SV_Capture","Adjust_Read_Lengths")))     x <- paste(x, "effective", sep = "_")
    x
}
svCapture_loadAssembly <- function(assemblyOptions, settings, rdsFile, ...){
    assemblyCache$get(
        'assembly', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            settings = settings$all(),
            rdsFile = rdsFile,
            internalUseSampleColumns = assemblyOptions$internalUseSampleColumns
        ), 
        createFn = function(...) {
            startSpinner(session, message = "loadAssembly.")
            assembly <- readRDS(rdsFile)
            startSpinner(session, message = "loadAssembly..")
            assembly$svs[, edgeType := svx_jxnType_altCodeToX(JXN_TYPE, "code")]
            assembly <- standardizeAssemblyColumns(assemblyOptions, assembly)
            startSpinner(session, message = "loadAssembly...")
            denominatorColumn <- svCapture_getDenominatorColumn(settings)
            assembly$samples$coverage <- as.integer(assembly$samples[[denominatorColumn]] / assembly$samples$sumTargetLens * 1000)
            assembly
        }
    )$value
}
svCapture_getGroups <- function(assemblyOptions, groupedProjectSamples, groupingCols, input, ...){
    assemblyCache$get(
        'groups', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            groupedProjectSamples,
            groupingCols,
            input$Data_Types
        ), 
        createFn = function(...) {
            startSpinner(session, message = "getGroups.")
            gps <- groupedProjectSamples
            gps$nSvs <- 0
            gps$svFrequency <- 0
            for(svType in input$Data_Types) {
                nSvs <- as.integer(gps[[svType]])
                gps$nSvs <- gps$nSvs + nSvs
                gps$svFrequency <- gps$svFrequency + nSvs / gps$coverage
            }
            startSpinner(session, message = "getGroups..")
            groups <- merge(
                aggegrateGroupSampleValues(gps, groupingCols, "svFrequency"),
                gps[, .(
                    coverage = sum(coverage),
                    nSvs = sum(nSvs)
                ), by = groupingCols]
            )
            setorderv(groups, groupingCols, order = -1L)
            startSpinner(session, message = "getGroups...")
            setAssemblyGroupLabels(groups, groupingCols)
        }
    )$value
}

#----------------------------------------------------------------------
# SVs found in the working samples, to populate distribution plots
#----------------------------------------------------------------------
svCapture_matchingAssemblySvs <- function(
    input,
    assemblyOptions, 
    assembly, groupedProjectSamples, groupingCols
){
    assembly <- assembly()
    gps <- groupedProjectSamples() 
    groupingCols <- groupingCols()   
    req(assembly, gps, groupingCols, input$Data_Types)
    assemblyCache$get(
        'matchingSvs', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            assembly$svs,
            gps,
            groupingCols,
            input$Data_Types
        ), 
        createFn = function(...) {
            startSpinner(session, message = "matching SVs.")
            svs <- gps[, {
                assembly$svs[
                    PROJECT == project & 
                    grepl(paste0(",", sample, ","), SAMPLES) & 
                    edgeType %in% svx_jxnType_longNameToX(input$Data_Types, "code")
                ]
            }, by = c("project", "sample", groupingCols)]
            startSpinner(session, message = "matching SVs..")
            svs <- setAssemblyGroupLabels(svs, groupingCols)
            stopSpinner(session)
            svs
        }
    )$value
}
svCapture_regroupedSvs <- function(
    input,
    assemblyOptions,
    assembly, groupedProjectSamples, groupingCols,
    conditions, groupLabels
){
    svCapture_matchingAssemblySvs(
        input,
        assemblyOptions,
        assembly, groupedProjectSamples, groupingCols
    )[
        groupLabel %in% groupLabels
    ] %>% regroupToUserConditions(
        groupingCols, conditions, groupLabels
    )
}

#----------------------------------------------------------------------
# capture target regions
#----------------------------------------------------------------------
svCapture_assemblyTargets <- function(assembly){
    assembly <- assembly() 
    req(assembly, assembly$targets)
    targets <- copy(assembly$targets)
    targets[, chromI := {
        x <- sub("chr", "", chrom)
        x <- ifelse(x == "X", 99, x)
        x <- ifelse(x == "Y", 100, x)
        as.integer(x)
    }]
    targets[, trackLabel := paste(name, paste0(chrom, geneStrand), sep = "\n")]
    targets[order(chromI, paddedStart)]
}
svCapture_parseDataByTarget <- function(
    input,
    assemblyOptions,
    assembly, groupedProjectSamples, groupingCols,
    conditions, groupLabels,
    xCol, dataFn
){
    targets <- svCapture_assemblyTargets(assembly)
    keyedTargets <- copy(targets)
    setkey(keyedTargets, name)
    x <- dataFn(conditions, groupLabels, targets)
    x$dt[, ":="(
        x = x$dt[[xCol]],
        trackLabel = unlist(keyedTargets[TARGET_REGION, trackLabel])
    )]
    x$trackLabels <- targets$trackLabel
    if(x$countGroups) x$groupCounts <- x$dt[, .N, by = .(groupLabel)]
    x$trackCounts <- x$dt[, .N, by = .(trackLabel)]
    x 
}
svCapture_xlimByTarget <- function(assembly, type){
    reactive({
        targets <- svCapture_assemblyTargets(assembly)
        x <- lapply(1:nrow(targets), function(i) switch(
            type,
            endpoints = targets[i, c(paddedStart + 1, paddedEnd)],
            coverage  = targets[i, {
                pad <- (regionEnd - regionStart) / 5
                c(regionStart + 1 - pad, regionEnd + pad)
            }]
        ) / 1e6)
        names(x) <- targets$trackLabel
        x
    })
}
svCapture_vLinesByTarget <- function(assembly){
    reactive({
        targets <- svCapture_assemblyTargets(assembly)
        x <- lapply(1:nrow(targets), function(i) targets[i, c(geneStart, geneEnd) / 1e6])
        names(x) <- targets$trackLabel
        x
    })
}
svCapture_vShadeByTarget <- function(assembly){
    reactive({
        targets <- svCapture_assemblyTargets(assembly)
        x <- lapply(1:nrow(targets), function(i) targets[i, c(regionStart, regionEnd) / 1e6])
        names(x) <- targets$trackLabel
        x
    })
}

#----------------------------------------------------------------------
# svCapture assembly plots
#----------------------------------------------------------------------
# SV frequencies
svCapture_svFrequenciesServer <- function(...){
    assemblyBarplotServer(..., ylab = "SV Frequency", nSD = 2)
}
#----------------------------------------------------------------------
# junction microhomology profiles
svCapture_microhomologyServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    plot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        defaultSettingValues = list(
            Density_Plot = list(
                Min_X_Value = -20,
                Max_X_Value =  20
            )
        ),
        xlab = "Insert Size (bp)",
        eventPlural = "SVs",        
        defaultBinSize = 1,
        v = c(seq(-50, 50, 5), -1, -2),
        x0Line = TRUE,
        dataFn = function(conditions, groupLabels) {
            svs <- svCapture_matchingAssemblySvs(
                input,
                assemblyOptions, 
                assembly, groupedProjectSamples, groupingCols
            )[
                groupLabel %in% groupLabels &
                N_SPLITS > 0 # omit gap-only junctions
            ] %>% regroupToUserConditions(
                groupingCols, conditions, groupLabels
            )
            svs$dt[, x := -MICROHOM_LEN]
            svs
        }
    )
}
#----------------------------------------------------------------------
# SV size profiles
svCapture_svSizesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    plot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "SV Size (log10 bp)",
        eventPlural = "SVs",        
        defaultBinSize = 0.1,
        v = 0:10,
        x0Line = FALSE,
        dataFn = function(conditions, groupLabels) {
            svs <- svCapture_matchingAssemblySvs(
                input,
                assemblyOptions, 
                assembly, groupedProjectSamples, groupingCols
            )[
                groupLabel %in% groupLabels
            ] %>% regroupToUserConditions(
                groupingCols, conditions, groupLabels
            )
            svs$dt[, x := log10(SV_SIZE)]
            svs
        }
    )
}
#----------------------------------------------------------------------
# junction endpoint profiles
svCapture_endpointsServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    getSvEndpoints <- function(conditions, groupLabels, targets){
        x <- svCapture_regroupedSvs(
            input,
            assemblyOptions,
            assembly, groupedProjectSamples, groupingCols,
            conditions, groupLabels
        )
        x$countGroups <- TRUE
        x$dt <- x$dt[, 
            .(
                POS = c(POS_1, POS_2) / 1e6,
                count = 1
            ), 
            by = .(TARGET_REGION, groupLabel)
        ]
        x
    }
    endpointsPlot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "SV Endpoint Coordinate (Mbp)",
        eventPlural = "ends",
        insideWidth = 2.5, 
        insideHeightPerBlock = 0.75,
        trackCols = "trackLabel",
        trackSameXLim = FALSE,
        trackSameYLim = TRUE,
        xlim = svCapture_xlimByTarget(assembly, "endpoints"),
        defaultBinSize = 10000 / 1e6,
        v = svCapture_vLinesByTarget(assembly),
        vShade = svCapture_vShadeByTarget(assembly),
        vColor = "grey40",
        x0Line = FALSE,
        y0Line = TRUE,
        # extraSettings = endpointSettings,
        # aggFn = sum,
        # aggCol = "count", 
        dataFn = function(conditions, groupLabels) {
            svCapture_parseDataByTarget(
                input,
                assemblyOptions,
                assembly, groupedProjectSamples, groupingCols,
                conditions, groupLabels,
                xCol = "POS", 
                dataFn = getSvEndpoints
            )
        }
    )
}
#----------------------------------------------------------------------
# target region coverage profiles
svCapture_coverageServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    coverageSettings <- list(
        Coverage = list(
            Aggregate_Coverage = list(
                type = "selectInput",
                choices = c("by sample","all samples together"),
                value = "by sample"
            )
        )
    )
    getTargetCoverage <- function(conditions, groupLabels, targets){
        assembly <- assembly()
        bins <- assembly$binnedCoverage
        env <- assembly$env 
        gps <- groupedProjectSamples()
        aggType <- coveragePlot$plot$settings$get("Coverage","Aggregate_Coverage")
        req(bins, gps)
        startSpinner(session, message = "getting coverage")
        projectSamples <- gps[, paste(project, sample, sep = "::")]
        samples        <- gps[, sample]
        dt <- do.call(rbind, lapply(targets$name, function(targetName){
            target <- targets[name == targetName]
            bins <- bins[[targetName]][, projectSamples]
            if(aggType == "by sample") {
                colnames(bins) <- samples
            } else {
                bins <- matrix(apply(bins, 1, sum, na.rm = TRUE), ncol = 1) 
                colnames(bins) <- "all samples"
            }
            do.call(rbind, lapply(1:ncol(bins), function(j){
                data.table(
                    TARGET_REGION = targetName,
                    groupLabel = colnames(bins)[j],
                    POS = (target$paddedStart + (1:nrow(bins) - 1) * env$COVERAGE_BIN_SIZE) / 1e6,
                    count = bins[, j]
                )
            }))
        }))
        if(aggType == "by sample") {
            groupLabels <- samples
            groupNs     <- 1
        } else {
            groupLabels <- "all samples"
            groupNs     <- length(samples)
        }
        stopSpinner(session)
        list(
            titleSuffix = NULL,
            groupLabels = groupLabels,
            groupCounts = data.table(groupLabel = groupLabels, N = groupNs),
            dt = dt,
            countGroups = FALSE
        )
    }
    coveragePlot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "SV Endpoint Coordinate (Mbp)",
        eventPlural = "samples",
        insideWidth = 2.5, 
        insideHeightPerBlock = 0.75,
        trackCols = "trackLabel",
        trackSameXLim = FALSE,
        trackSameYLim = TRUE,
        xlim = svCapture_xlimByTarget(assembly, "coverage"),
        defaultBinSize = 1000 / 1e6,
        v = svCapture_vLinesByTarget(assembly),
        vShade = svCapture_vShadeByTarget(assembly),
        vColor = "grey40",
        x0Line = FALSE,
        y0Line = TRUE,
        extraSettings = coverageSettings,
        aggFn = sum,
        aggCol = "count", 
        dataFn = function(conditions, groupLabels) {
            svCapture_parseDataByTarget(
                input,
                assemblyOptions,
                assembly, groupedProjectSamples, groupingCols,
                conditions, groupLabels,
                xCol = "POS", 
                dataFn = getTargetCoverage
            )
        }
    )
}
#----------------------------------------------------------------------
# junction base usage profiles
svCapture_junctionBasesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    # junctionBasesSettings <- list(
    #     Junction_Bases = list(
    #         Aggregate_Coverage = list(
    #             type = "selectInput",
    #             choices = c("by sample","all samples together"),
    #             value = "by sample"
    #         )
    #     )
    # )
    baseValues <- list(A = 1L, C = 2L, G = 3L, T = 4L)
    baseLookup <- names(baseValues)
    getBaseUsageProfile <- function(conditions, groupLabels){
        assembly_ <- assembly()
        bases <- assembly_$baseUsageProfile # [SV_side, position], value = numeric base, rownames = paste(project, svId, side)
        env <- assembly_$env 
        svs <- svCapture_regroupedSvs(
            input,
            assemblyOptions,
            assembly, groupedProjectSamples, groupingCols,
            conditions, groupLabels
        )$dt
        req(bases, svs)
        svSvKeys <- svs[, paste(project, SV_ID, sep = "::")]
        basesSvKeys <- sapply(rownames(bases), function(x) paste(strsplit(x, "::")[[1]][1:2], collapse = "::"))
        I <- which(basesSvKeys %in% svSvKeys)
        bases <- bases[I, ]
        list(
            groupLabels = baseValues,
            groupCounds = {
                x <- as.vector(bases)
                x <- aggregate(x, list(x), length)
                names(x) <- c("groupLabel", "N")
                x
            },
            dt = do.call(rbind, lapply(1:ncol(bases), function(pos){
                data.table(
                    x = pos - env$BASE_USAGE_SPAN - 1,
                    groupLabel = baseLookup[bases[, pos]]
                )
            }))
        )
    # junctionBaseUsageProfile = junctionBaseUsageProfile,
    # targetBaseUsageProfile = targetBaseUsageProfile,
    }
    baseUsagePlot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "Position Relative to SV Junction (bp)",
        eventPlural = "positions",
        insideWidth = 2.5, 
        insideHeightPerBlock = 1.5,
        xlim = reactive({
            BASE_USAGE_SPAN <- assembly()$env$BASE_USAGE_SPAN
            c(-BASE_USAGE_SPAN, BASE_USAGE_SPAN)
        }),
        defaultBinSize = 1,
        x0Line = TRUE,
        y0Line = TRUE,
        # extraSettings = junctionBasesSettings,
        dataFn = getBaseUsageProfile
    )
}
#----------------------------------------------------------------------
svx_getInsertionsTitle <- function(assembly, settings, nSvs, family = "Plot"){
    title <- settings$get("Plot", "Title", NA)
    nSvs <- paste(commify(nSvs), "SVs")
    if(!isTruthy(title) && !is.null(assembly)) title <- assembly$env$DATA_NAME
    if(isTruthy(title)) paste0(title, " (", nSvs, ")")
    else nSvs
}
svx_getInsertionsProperty <- function(assembly, settings, assemblyKey, settingsKey, default = NA){
    if(is.null(assembly) || is.null(assembly$env[[assemblyKey]])){
        settingsFamily <- if("Junction_Properties" %in% names(settings)) "Junction_Properties" else "Insertions"
        if(is.null(settings[[settingsFamily]])) return(default)
        settings[[settingsFamily]]()[[settingsKey]]$value
    } else {
        assembly$env[[assemblyKey]]
    }
}
svx_getInsertionPlotMetadata <- function(assembly, settings){
    minInsertionSize <- svx_getInsertionsProperty(assembly, settings, "MIN_INSERTION_SIZE", "Min_Insertion_Size", default = 1)
    maxInsertionSize <- svx_getInsertionsProperty(assembly, settings, "MAX_INSERTION_SIZE", "Max_Insertion_Size", default = 20)
    searchSpace      <- svx_getInsertionsProperty(assembly, settings, "INSERTION_SEARCH_SPACE", "Insertion_Search_Space", default = 500)
    minTemplateSize  <- svx_getInsertionsProperty(assembly, settings, "MIN_TEMPLATE_SIZE", "Min_Template_Size",   default = 7)
    foldback     <- CONSTANTS$plotlyColors$blue
    palindrome   <- CONSTANTS$plotlyColors$red
    crossJxn     <- CONSTANTS$plotlyColors$green
    slippage     <- CONSTANTS$plotlyColors$orange
    strandSwitch <- CONSTANTS$plotlyColors$purple
    other        <- CONSTANTS$plotlyColors$black
    typeColors <- list(
        foldback     = foldback, 
        palindrome   = palindrome,
        crossJxn     = crossJxn,
        slippage     = slippage,
        strandSwitch = strandSwitch, 
        other        = other
    )
    typeLabels <- c(
        foldback     = "foldback", 
        palindrome   = "palindrome",
        crossJxn     = "cross-junction",
        slippage     = "slippage",
        strandSwitch = "strand-switch", 
        other        = "other"
    )
    typeStrands <- c(
        foldback     = -1, 
        palindrome   = -1,
        crossJxn     =  1,
        slippage     =  1,
        strandSwitch = -1, 
        other        =  1
    )
    list(
        minInsertionSize = minInsertionSize,
        maxInsertionSize = maxInsertionSize,
        minTemplateSize  = minTemplateSize,
        searchSpace      = searchSpace,
        totalSearchSpace = searchSpace * 2 * 2 * 2, # two breakpoints searched on both sides of the breakpoint on two strands 
        pValueThreshold  = svx_getInsertionsProperty(NULL, settings, NULL, "P_Value_Threshold", default = 0.05),
        maxDist          = svx_getInsertionsProperty(NULL, settings, NULL, "Max_Plotted_Distance",  default = 2),
        vSpacing         = svx_getInsertionsProperty(NULL, settings, NULL, "Distance_Grid_Spacing", default = 5),
        nCombinations = sapply(1:maxInsertionSize, function(insertionSize){
            flankUHom <- if(insertionSize >= minTemplateSize) 1L
                        else as.integer(ceiling((minTemplateSize - insertionSize) / 2))
            4 ** (insertionSize + 2 * flankUHom)
        }),
        typeColors  = typeColors,
        typeLabels  = typeLabels,
        typeStrands = typeStrands
    )
}
svx_insertions_addTopLegend <- function(plot, md){
    plot$addLegend(
        legend = md$typeLabels[names(md$typeLabels) != "other"],
        col =    unlist(md$typeColors[names(md$typeColors) != "other"]),
        x = "top",
        horiz = TRUE,
        x.intersp = 0,
        pch = 19,
        pt.cex = 1 
    )
}
#----------------------------------------------------------------------
svx_plotInsertions_yield <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    yield <- svs[, .(
        nCombinations = md$nCombinations[-MICROHOM_LEN], # number of possible unique template sequences for an insertion size 
        nSvs = .N,                                    # number of trials ...
        nFound = sum(templateType != "notFound")      # number of successes ...
    ), keyby = MICROHOM_LEN][,                        # ... by insert size
        mu := md$totalSearchSpace * (1 / nCombinations)       # Poisson mean, i.e., expected target density ...
    ][, ":="(
        successRate = nFound / nSvs, # fraction of SVs whose search sequence was found
        trialSuccessProb = 1 - dpois(0, mu) # probably of finding at least one match
    )][, 
        pValue := 1 - pbinom(nFound - 1, nSvs, trialSuccessProb)
    ]
    plot$initializeFrame(
        xlim = range(-yield$MICROHOM_LEN),
        ylim = c(0, min(100, max(yield$successRate, yield$trialSuccessProb) * 110)),
        xlab = "Insertion Size (bp)",
        ylab = "% Found Templates",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(yield$nSvs))
    )  
    plot$addLines(
        x = -yield$MICROHOM_LEN,
        y = yield$trialSuccessProb * 100,
        col = CONSTANTS$plotlyColors$blue
    )     
    plot$addPoints(
        x = -yield$MICROHOM_LEN,
        y = yield$successRate * 100,
        col = ifelse(yield$pValue < md$pValueThreshold, 
                        CONSTANTS$plotlyColors$red, 
                        CONSTANTS$plotlyColors$black)
    )
    plot$addLegend(
        x = "topright",
        legend = c(
            paste("p", "<",  md$pValueThreshold),
            paste("p", ">=", md$pValueThreshold),
            "expected"
        ),
        col = c(
            CONSTANTS$plotlyColors$red, 
            CONSTANTS$plotlyColors$black,
            CONSTANTS$plotlyColors$blue
        ),
        pch = c(19, 19, NA),
        pt.cex = c(1, 1, NA),
        lty = c(NA, NA, 1),
        lwd = c(NA, NA, 1)
    )
}
svx_plotInsertions_map <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    svs <- svs[templateType != "notFound"]
    svs[, ":="(
        x = switch(
            templateBreakpointN,
            switch(
                templateType,
                foldback      = -templateDistance,
                palindrome    = -templateDistance,
                crossJxn      = -templateDistance,
                slippage      =  templateDistance,
                strandSwitch  =  templateDistance,
                other         =  templateDistance
            ),
            switch(
                templateType,
                foldback      =  templateDistance,
                palindrome    =  templateDistance,
                crossJxn      =  templateDistance,
                slippage      = -templateDistance,
                strandSwitch  = -templateDistance,
                other         = -templateDistance
            )
        ),
        y = switch(
            templateBreakpointN,
            if(templateIsRc) 3 else 4,
            if(templateIsRc) 1 else 2
        )
    ), by = c("PROJECT", "SV_ID")]
    plot$initializeFrame(
        xlim = c(-md$maxDist, md$maxDist),
        ylim = c(0.5, 6),
        xlab = "Insertion Template Distance from Junction (bp)",
        ylab = "",
        yaxt = "n",
        bty = "n",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(nrow(svs)))
    )

    lastLine <- floor(md$maxDist / md$vSpacing) * md$vSpacing
    for(x in seq(-lastLine, lastLine, md$vSpacing)){
        lines(c(x, x), c(0, 4.75), col = CONSTANTS$plotlyColors$grey)
    }
    boldLwd <- 2.5
    for(y in 3:4){ # left side of junction, thick line drawn on top two lines
        lines(c(-lastLine, 0), c(y, y), col = CONSTANTS$plotlyColors$black, lwd = boldLwd)  
        lines(c(0, lastLine),  c(y, y), col = CONSTANTS$plotlyColors$grey)   
    }
    for(y in 1:2){ # right side of junction, thick line drawn on bottom two lines
        lines(c(-lastLine, 0), c(y, y), col = CONSTANTS$plotlyColors$grey)  
        lines(c(0, lastLine),  c(y, y), col = CONSTANTS$plotlyColors$black, lwd = boldLwd)   
    }
    lines(c(0, 0), c(1, 4), col = CONSTANTS$plotlyColors$black, lwd = boldLwd)
    
    text(-md$maxDist, 5.25, "Left Side of SV",  pos = 4)
    text( md$maxDist, 5.25, "Right Side of SV", pos = 2)

    plot$addPoints(
        x = jitter(svs$x, amount = 0.25), # jitter to make all points visible
        y = jitter(svs$y, amount = 0.25),
        col = sapply(svs$templateType, function(tt) md$typeColors[[tt]])
    )
    svx_insertions_addTopLegend(plot, md)
}
svx_plotInsertions_histogram <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    svs <- svs[templateType != "notFound"]
    d <- svs[, .(
        templateType,
        templateBreakpointN,
        pos = templateStartPos:templateEndPos - md$searchSpace - if(templateBreakpointN == 1) 0 else 1
    ), by = .(PROJECT, SV_ID)] %>% dcast(
        templateBreakpointN + pos ~ templateType, 
        fun.aggregate = length, 
        fill = 0
    )
    posRange <- min(d$pos):max(d$pos)
    d <- Reduce(
        function(dt1, dt2) merge(dt1, dt2, by = c("templateBreakpointN", "pos"), all.x = TRUE, all.y = TRUE),
        list(
            data.table(templateBreakpointN = 1, pos = posRange), 
            data.table(templateBreakpointN = 2, pos = posRange),
            d
        )
    )
    d[is.na(d)] <- 0
    templateTypes <- intersect(names(md$typeColors), names(d))
    dstr(templateTypes)
    dstr(d)
    plot$initializeFrame(
        xlim = c(-md$maxDist, md$maxDist),
        ylim = c(0, d[, max(.SD, na.rm = TRUE) * 1.5, .SDcols = templateTypes]),
        xlab = "Distance from Junction (bp)",
        ylab = "# of Templates",
        # yaxt = "n",
        bty = "n",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(nrow(svs)))
    )
    # abline(v = seq(vSpacing, maxDist - 1, vSpacing), col = CONSTANTS$plotlyColors$grey)
    abline(v = 0, col = CONSTANTS$plotlyColors$grey)
    lwd <- 1.5 # plot$settings$Points_and_Lines()$Line_Width$value
    for(bpn in 1:2) for(type in templateTypes){
        dd <- d[templateBreakpointN == bpn]
        lines(dd$pos, dd[[type]], col = md$typeColors[[type]], lwd = lwd, lty = bpn)
    }
    svx_insertions_addTopLegend(plot, md)
}
svx_plotInsertions_pileup <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    svs <- svs[templateType != "notFound"]
    templateTypes <- names(md$typeColors)
    svs <- svs[
        templateType %in% templateTypes, 
        .(
            templateBreakpointN,
            templateType,
            templateStartFm,
            templateEndFm,
            series = paste(templateBreakpointN, templateType, sep = ":"),
            # series = if(templateType == "other") 0L else templateBreakpointN,
            startPos = templateStartPos - md$searchSpace - if(templateBreakpointN == 1) 0L else 1L,
            endPos   = templateEndPos   - md$searchSpace - if(templateBreakpointN == 1) 0L else 1L,
            size     = templateEndPos - templateStartPos + 1
        ), 
        by = .(PROJECT, SV_ID)
    ][, ":="(
        insStartPos = startPos + templateStartFm,
        insEndPos   = endPos   - templateEndFm
    )][
        abs(startPos) <= md$maxDist | 
        abs(endPos)   <= md$maxDist
    ][
        order(series, -size)
    ][,
        plotRow := 1:.N,
        by = .(series)
    ]
    maxNSvs <- svs[, .N, by = .(series)][, max(N)]
    plot$initializeFrame(
        xlim = c(-md$maxDist, md$maxDist),
        ylim = c(0, maxNSvs * 4 * 1.1),
        xlab = "Distance from Junction (bp)",
        ylab = "", # "# of Templates",
        yaxt = "n",
        bty = "n",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(nrow(svs)))
    )
    abline(v = 0, col = CONSTANTS$plotlyColors$grey)
    abline(h = c(1,3) * maxNSvs, col = CONSTANTS$plotlyColors$grey)
    text( md$maxDist / 2, 3.5 * maxNSvs, "Breakpoint #1")
    text(-md$maxDist / 2, 0.5 * maxNSvs, "Breakpoint #2")
    lwd <- 0.75 # plot$settings$Points_and_Lines()$Line_Width$value
    svs[, {
        yStrand <- md$typeStrands[[templateType[1]]]
        chunkScalar <- if(templateBreakpointN[1] == 1) 3 else 1
        yOffset <- chunkScalar * maxNSvs
        for(svI in seq_len(.N)){
            .SD[svI, {
                y <- plotRow * yStrand + yOffset
                lines(c(startPos, endPos), c(y, y), col = md$typeColors[[templateType]],  lwd = lwd)
            }]
        }
    }, by = .(series)]
    svx_insertions_addTopLegend(plot, md)
}
svx_plotInsertions_microhomologyLengths <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    svs <- svs[templateType != "notFound"]
    xylim <- c(0.5, 8.5)
    plot$initializeFrame(
        xlim = xylim,
        ylim = xylim,
        xlab = "Lef microhomology length (bp)",
        ylab = "Right microhomology length (bp)",
        title = svx_getInsertionsTitle(assembly, plot$settings, nrow(svs))
    )  
    levelPlotSettings <- list(
        Level_Plot = list(
            Max_Z_Value = list(
                value = "auto"
            ),
            Level_Palette = list(
                value = "seq Blues"
            ),
            Legend_Digits = list(
                value = 1
            )
        )
    )
    mdiLevelPlot(
        svs[, .(
            x = templateStartFm,
            y = templateEndFm,
            templateType
        )],
        xlim = xylim,
        xinc = 1,
        ylim = xylim,
        yinc = 1,
        z.fn = length,   # function applied to z.columnumn, per grid spot, to generate the output color
        z.column = "templateType", # the column in dt passed to z.fn, per grid spot
        settings = levelPlotSettings, # a settings object from the enclosing staticPlotBox, or any list compatible with mdiLevelPlotSettings
        legendTitle = "# of SVs"
    )
}
svx_plotInsertions_nTemplateInstances <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    # svs <- svs[templateType != "notFound"]
    svs[is.na(templateInstances), templateInstances := 0]
    density <- svs[, .(
        nSvs = .N
    ), keyby = .(templateInstances)]
    plot$initializeFrame(
        xlim = c(min(density$templateInstances) - 0.5, max(density$templateInstances) + 0.5),
        ylim = c(0, max(density$nSvs) * 1.05),
        xlab = "# of Found Templates",
        ylab = "# of SVs",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(density$nSvs))
    )     
    plot$addPoints(
        x = density$templateInstances,
        y = density$nSvs,
        col = CONSTANTS$plotlyColors$blue,
        typ = "h",
        lwd = 6
    )
}
# svx_plotInsertions_templateDistances <- function(plot, svs, assembly = NULL){
#     md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
#     svs <- svs[templateType != "notFound"]
#     density <- svs[, .(
#         nSvs = .N
#     ), keyby = .(templateType, templateDistance)]
#     plot$initializeFrame(
#         xlim = c(0, md$maxDist),
#         ylim = c(0, max(density$nSvs) * 1.1),
#         xlab = "Template Distance from Junction",
#         ylab = "# of SVs",
#         title = svx_getInsertionsTitle(assembly, plot$settings, sum(density$nSvs))
#     )     
#     for(type in density$templateType){
#         dd <- density[templateType == type]
#         plot$addLines(
#             x = dd$templateDistance,
#             y = dd$nSvs,
#             col = md$typeColors[[type]]
#         )        
#     }
#     svx_insertions_addTopLegend(plot, md)
# }
svCapture_insertionTemplatesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    settings <- c(assemblyPlotFrameSettings, list(
        Insertions = list(
            Max_Plotted_Distance = list(
                type = "numericInput",
                value = 25
            ),
            Distance_Grid_Spacing = list(
                type = "numericInput",
                value = 5
            ),
            Plot_Insertions_As = list(
                type = "selectInput",
                choices = c(
                    "yield",
                    "map",
                    "histogram",
                    "pileup",
                    "microhomologyLengths",
                    "nTemplateInstances"
                    # ,
                    # "templateDistances"
                ),
                value = "pileup"
            ),
            P_Value_Threshold = list(
                type = "numericInput",
                value = 0.05
            ),
            Allow_Multiply_Matched_SVs = list(
                type = "checkboxInput",
                value = TRUE
            ),
            High_Probability_Only = list(
                type = "checkboxInput",
                value = FALSE
            )
        )
    ))
    mars <- list(
        yield = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        map = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$nullMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        histogram = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        pileup = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        microhomologyLengths = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar
        ),
        nTemplateInstances = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        templateDistances = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        )
    )
    dims <- list(
        yield = list(
            width  = 1.5, 
            height = 1.25
        ),
        map = list(
            width  = 4, 
            height = 2
        ),
        histogram = list(
            width  = 4, 
            height = 1.5
        ),
        pileup = list(
            width  = 4, 
            height = 3
        ),
        microhomologyLengths = list(
            width  = 1.5, 
            height = 1.5
        ),
        nTemplateInstances = list(
            width  = 1.5, 
            height = 1.5
        ),
        templateDistances = list(
            width  = 1.5, 
            height = 1.5
        )
    )
    getDim <- function(plotAs, key, option){
        dim <- trimws(assemblyPlot$plot$settings$get("Plot",option))
        if(!isTruthy(dim) || dim == "" || dim == "auto") dims[[plotAs]][[key]]
        else dim
    }
    assemblyPlot <- assemblyPlotBoxServer( 
        id, session, input, output, 
        isProcessingData,
        groupingCols, groups,
        dataFn = function(conditions, groupLabels) {
            svCapture_matchingAssemblySvs(
                input,
                assemblyOptions, 
                assembly, groupedProjectSamples, groupingCols
            )[
                groupLabel %in% groupLabels & 
                !is.na(templateType) # thus, all insertion junctions from selected samples that were analyzed (even if template not found)
            ]
        },  
        plotFrameFn = function(data) {
            plotAs <- assemblyPlot$plot$settings$get("Insertions","Plot_Insertions_As")
            list(
                frame = getAssemblyPlotFrame(
                    plot = assemblyPlot$plot, 
                    insideWidth  = getDim(plotAs, "width",  "Width_Inches"), 
                    insideHeight = getDim(plotAs, "height", "Height_Inches"), 
                    mar = mars[[plotAs]]
                ),
                mar = mars[[plotAs]]
            )
        },
        plotFn = function(plotId, dataReactive, plotFrameReactive) staticPlotBoxServer(
            plotId,
            settings = settings, 
            size = "m",
            Plot_Frame = reactive({ plotFrameReactive()$frame }),
            create = function() {
                d <- dataReactive()
                plotAs <- assemblyPlot$plot$settings$get("Insertions","Plot_Insertions_As")
                assembly <- assembly()
                req(d, d$data, plotAs, assembly)
                startSpinner(session, message = paste("rendering", id))
                par(mar = plotFrameReactive()$mar)
                svs <- d$data
                allowMultiple <- assemblyPlot$plot$settings$get("Insertions","Allow_Multiply_Matched_SVs")
                if(!allowMultiple) svs <- svs[
                    is.na(templateInstances) | 
                    templateInstances <= 1
                ]
                highProbability <- assemblyPlot$plot$settings$get("Insertions","High_Probability_Only")
                if(highProbability) svs <- svs[
                    !is.na(templateType) &
                    templateType != "other" & 
                    templateDistance <= 20
                ]
                get(paste("svx_plotInsertions", plotAs, sep = "_"))(assemblyPlot$plot, svs, assembly)
                stopSpinner(session)
            }
        )
    )
}


#   [1] "ACTION_DIR"               "GENOMEX_MODULES_DIR"     
# mdi-web-server-app-server-1              |  [3] "INPUT_DIR"                "SAMPLES_TABLE"           
# mdi-web-server-app-server-1              |  [5] "OUTPUT_DIR"               "DATA_NAME"               
# mdi-web-server-app-server-1              |  [7] "DATA_FILE_PREFIX"         "ASSEMBLE_PREFIX"         
# mdi-web-server-app-server-1              |  [9] "N_CPU"                    "MIN_FAMILY_SIZE"         
# mdi-web-server-app-server-1              | [11] "MIN_SV_SIZE"              "MAX_SV_SIZE"             
# mdi-web-server-app-server-1              | [13] "MIN_SAMPLES"              "MAX_SAMPLES"             
# mdi-web-server-app-server-1              | [15] "MIN_MOLECULES"            "MAX_MOLECULES"           
# mdi-web-server-app-server-1              | [17] "MIN_SPLIT_READS"          "COVERAGE_BIN_SIZE"       
# mdi-web-server-app-server-1              | [19] "MIN_INSERTION_SIZE"       "MAX_INSERTION_SIZE"      
# mdi-web-server-app-server-1              | [21] "MIN_TEMPLATE_SIZE"        "INSERTION_SEARCH_SPACE"  
# mdi-web-server-app-server-1              | [23] "BASE_USAGE_SPAN"          "MAX_SHARED_PROPER"       
# mdi-web-server-app-server-1              | [25] "WARN_ONLY"                "INCLUDE_CLIPS"           
# mdi-web-server-app-server-1              | [27] "FORCE_COVERAGE"           "FIND_INSERTION_TEMPLATES"
# mdi-web-server-app-server-1              | [29] "PROFILE_BASE_USAGE"  