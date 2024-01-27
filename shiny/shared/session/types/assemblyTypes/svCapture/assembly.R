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
            valueColumn <- "svFrequency"
            for(svType in input$Data_Types) {
                nSvs <- as.integer(gps[[svType]])
                gps$nSvs <- gps$nSvs + nSvs
                gps$svFrequency <- gps$svFrequency + nSvs / gps$coverage
                dtvc <- paste(valueColumn, svType, sep = "__")
                gps[[dtvc]] <- nSvs / gps$coverage
            }
            startSpinner(session, message = "getGroups..")
            isSingleGroup <- length(groupingCols) == 0
            bBy <- if(isSingleGroup) NULL else groupingCols
            a <- aggegrateGroupSampleValues(gps, groupingCols, valueColumn, input)
            b <- gps[, .(
                coverage = sum(coverage),
                nSvs = sum(nSvs)
            ), by = bBy]
            groups <- if(isSingleGroup) cbind(a, b) else merge(a, b)
            if(!isSingleGroup) setorderv(groups, groupingCols, order = -1L)
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
    req(assembly, gps, input$Data_Types)
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
    ) # return object with key and value to support assemblyCache$get cascading
}
svCapture_invertInsertionsByGene <- function(svs, assembly, assemblyOptions, targets){
    assembly <- assembly()
    req(svs, assembly, assembly$env$INSERTION_SEARCH_SPACE)
    assemblyCache$get(
        'invertInsertions', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            svs$key,
            targets,
            assembly$env$INSERTION_SEARCH_SPACE
        ), 
        createFn = function(...) {
            startSpinner(session, message = "inverting insertions")
            ss <- 2 * assembly$env$INSERTION_SEARCH_SPACE
            targets <- copy(targets)
            setkey(targets, name)
            svs <- svs$value
            svs[!is.na(templateType), ":="(
                geneStrand = targets[TARGET_REGION, geneStrand],
                templateIsInverted = FALSE
            )]
            svs[!is.na(geneStrand) & geneStrand == "-", ":="( # top strand AND any SVs where TARGET_REGION not in target$name pass as is
                templateStartPos = ss - templateEndPos + 1, # does this need to use searchSpace? if so, need assembly()
                templateStartFm  = templateEndFm,
                templateStartIsRetained = templateEndIsRetained,
                templateEndPos   = ss - templateStartPos + 1,
                templateEndFm    = templateStartFm,
                templateEndIsRetained = templateStartIsRetained,
                templateBreakpointN = 3 - templateBreakpointN,
                templateIsInverted = TRUE
            )]
            stopSpinner(session)
            svs
        }
    )$value
}
svCapture_invertArcsByGene <- function(svs, assembly, assemblyOptions, targets){
    assembly <- assembly()
    req(svs, assembly)
    assemblyCache$get(
        'invertArcs', 
        permanent = TRUE,
        from = "ram",

        #######################
        # create = "once",  #assemblyOptions$cacheCreateLevel,
        create = assemblyOptions$cacheCreateLevel,

        keyObject = list(
            svs$key,
            targets
        ), 
        createFn = function(...) {
            startSpinner(session, message = "inverting arcs")
            svs <- svs$value
            targets <- copy(targets)
            setkey(targets, name)
            startSpinner(session, message = "inverting intraTarget")
            intraTargetSvs <- svs[TARGET_CLASS %in% c("TT","TA") & TARGET_REGION %in% targets$name]
            intraTarget <- if(nrow(intraTargetSvs) == 0) NULL else intraTargetSvs[, # always Del, Dup or Inv
                {
                    geneCenter <- targets[TARGET_REGION, geneStart + (geneEnd - geneStart) / 2]
                    invert     <- targets[TARGET_REGION, geneStrand == "-"]
                    pos <- sort(c(POS_1, POS_2))
                    offsets <- pos - geneCenter
                    sides <- c(SIDE_1, SIDE_2)[order(c(POS_1, POS_2))]
                    .(
                        TARGET_CLASS  = TARGET_CLASS,
                        TARGET_REGION = TARGET_REGION,
                        JXN_TYPE      = JXN_TYPE,
                        OFFSET_1 = if(invert) -offsets[2] else offsets[1], # for intraTarget, OFFSET_1 is gene-proximal
                        OFFSET_2 = if(invert) -offsets[1] else offsets[2],
                        SIDE_1   = if(invert) { if(sides[2] == "L") "R" else "L" } else sides[1],
                        SIDE_2   = if(invert) { if(sides[1] == "L") "R" else "L" } else sides[2],
                        groupLabel = groupLabel
                    )
                }, 
                by = .(PROJECT, SV_ID)
            ]
            startSpinner(session, message = "inverting interTarget")
            interTargetSvs <- svs[TARGET_CLASS %in% c("tt","ta") & grepl(",", TARGET_REGION)]
            interTarget <- if(nrow(interTargetSvs) == 0) NULL else interTargetSvs[, # often Trans, but could be Del, Dup, Inv if two targets are on the same chromosome
                {
                    targets_ <- rbind(
                        targets[
                            CHROM_1 == chrom & 
                            POS_1 >= paddedStart & 
                            POS_1 <= paddedEnd
                        ][1],
                        targets[
                            CHROM_2 == chrom & 
                            POS_2 >= paddedStart & 
                            POS_2 <= paddedEnd
                        ][1]
                    )[ , i := 1:2]
                    i <- targets_[order(chromI, geneStart), i]
                    geneCenters <- targets_[i, geneStart + (geneEnd - geneStart) / 2]
                    inverts     <- targets_[i, geneStrand == "-"]
                    offsets <- c(POS_1, POS_2)[i] - geneCenters
                    sides <- c(SIDE_1, SIDE_2)[i]
                    .(
                        TARGET_CLASS  = TARGET_CLASS,
                        TARGET_REGION = TARGET_REGION,
                        JXN_TYPE      = JXN_TYPE,
                        OFFSET_1 = if(inverts[1]) -offsets[1] else offsets[1], # for translocation, OFFSET_1 is for the lower sorting gene
                        OFFSET_2 = if(inverts[2]) -offsets[2] else offsets[2],
                        SIDE_1   = if(inverts[1]) { if(sides[1] == "L") "R" else "L" } else sides[1],
                        SIDE_2   = if(inverts[2]) { if(sides[2] == "L") "R" else "L" } else sides[2],
                        groupLabel = groupLabel
                    )
                }, 
                by = .(PROJECT, SV_ID)
            ]
            stopSpinner(session)
            list(
                intraTarget = intraTarget,
                interTarget = interTarget,
                maxOffset = targets[, max(paddedEnd - paddedStart) / 2]
            )
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
    )$value[
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
svCapture_xlimByTarget <- function(assembly, type, assemblyPlot){
    reactive({
        targets <- svCapture_assemblyTargets(assembly)
        x <- lapply(1:nrow(targets), function(i) {
            paddedRegion <- targets[i, c(paddedStart + 1, paddedEnd)]
            targetRegion <- targets[i, {
                pad <- (regionEnd - regionStart) / 5
                c(regionStart + 1 - pad, regionEnd + pad)
            }]
            x <- switch(
                type,
                endpoints = {
                    normalize <- assemblyPlot$plot$settings$get("Endpoints","Normalize_to_Coverage")
                    if(normalize) targetRegion else paddedRegion
                },
                coverage  = targetRegion
            ) / 1e6
            if(targets[i, geneStrand == "-"]) x <- rev(x) # thus, all target-track plots are drawn in transcribed order, left to right
            x
        })
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
# ----------------------------------------------------------------------
svCapture_svFrequenciesServer <- function(...){
    assemblyBarplotServer(..., ylab = "SV Frequency", nSD = 2)
}
#----------------------------------------------------------------------
# junction microhomology profiles
# ----------------------------------------------------------------------
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
            )$value[
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
# ----------------------------------------------------------------------
splitSvsByJxnType <- function(svs, plot, family, groupingCols){
    if(!isTruthy(plot$settings$get(family,"Split_by_SV_Type"))) return(svs)
    groupingCols <- groupingCols()
    svTypes <- svs$dt[, paste("SV_type =", tolower(svx_jxnType_altCodeToX(JXN_TYPE, "longName")))]
    svs$dt[, groupLabel := if(length(groupingCols) == 0) svTypes else paste(groupLabel, svTypes, sep = " | ")]
    svs$groupLabels <- unique(svs$dt$groupLabel)
    svs$groupCounts <- svs$dt[, .N, by = .(groupLabel)]
    svs
}
svCapture_svSizesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    sizeSettings <- list(
        SV_Sizes = list(
            Split_by_SV_Type = list(
                type = "checkboxInput",
                value = FALSE
            )
        )
    )
    sizesPlot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "SV Size (log10 bp)",
        eventPlural = "SVs",
        extraSettings = sizeSettings,
        defaultBinSize = 0.1,
        v = 1:10,
        vColor = "grey80",
        insideHeightPerBlock = 0.75,
        groupV = function(data) sapply(data$groupLabels, function(groupLabel_){
            data$dt[groupLabel == groupLabel_, median(x, na.rm = TRUE)]
        }, simplify = FALSE, USE.NAMES = TRUE),
        x0Line = FALSE,
        dataFn = function(conditions, groupLabels) {
            svs <- svCapture_matchingAssemblySvs(
                input,
                assemblyOptions, 
                assembly, groupedProjectSamples, groupingCols
            )$value[
                groupLabel %in% groupLabels
            ] %>% regroupToUserConditions(
                groupingCols, conditions, groupLabels
            ) %>% splitSvsByJxnType(sizesPlot$plot, "SV_Sizes", groupingCols)
            svs$dt[, x := log10(SV_SIZE)]
            svs
        }
    )
}
#----------------------------------------------------------------------
# junction endpoint profiles
# ----------------------------------------------------------------------
svCapture_getTargetCoverages <- function(
    assembly, groupedProjectSamples, groupingCols, aggType, plot, 
    conditions, groupLabels
){
    targets <- svCapture_assemblyTargets(assembly)
    assembly <- assembly()
    binCov <- assembly$binnedCoverage # named list, key = regionKey, value = int[bin, sample], colnames = project::sample
    binSize <- assembly$env$COVERAGE_BIN_SIZE
    req(binCov)
    startSpinner(session, message = "binning coverage")
    groupedProjectSamples <- groupedProjectSamples()
    groupingCols <- groupingCols()
    gps <- setAssemblyGroupLabels(groupedProjectSamples, groupingCols) 
    req(gps)
    gps[, projectSample := paste(project, sample, sep = "::")]
    projectSamples <- gps[groupLabel %in% groupLabels, projectSample]
    if(is.null(aggType)) aggType <- plot$settings$get("Coverage","Aggregate_Coverage")
    mergeGroups <- aggType == "all samples together"
    dt <- do.call(rbind, lapply(targets$regionKey, function(regionKey_){
        target <- targets[regionKey == regionKey_]
        ps <- projectSamples[projectSamples %in% colnames(binCov[[regionKey_]])]
        binCov <- binCov[[regionKey_]][, ps, drop = FALSE]
        if(mergeGroups) {
            binCov <- matrix(apply(binCov, 1, sum, na.rm = TRUE), ncol = 1) 
            colnames(binCov) <- "all samples"
        }
        x <- as.integer(target$paddedStart + (1:nrow(binCov) - 1) * binSize) / 1e6
        do.call(rbind, lapply(colnames(binCov), function(col){
            data.table(
                regionKey = target$regionKey,
                track = target$trackLabel,
                group = if(mergeGroups) col else gps[projectSample == col, groupLabel],
                x = x,
                y = binCov[, col] / 1e3
            )
        }))
    }))[, .(y = sum(y, na.rm = TRUE)), .(track, group, x)]
    if(aggType == "by group") {
        groupLabels <- groupLabels
        groupNs     <- sapply(groupLabels, function(x) gps[groupLabel == x, .N])
    } else {
        groupLabels <- "all samples"
        groupNs     <- length(projectSamples)
    }
    list( 
        titleSuffix = NULL,
        groupLabels = groupLabels,
        groupCounts = data.table(groupLabel = groupLabels, N = groupNs), 
        countGroups = FALSE,
        dt = list( # this dt level "tricks" assemblyDensityPlotServer
            X_Bin_Size = binSize / 1e6, 
            totalN = length(projectSamples),
            trackLabels = targets$trackLabel,
            dt = dt
        ),
        trackLabels = targets$trackLabel
    )
}
svCapture_endpointsServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    endpointSettings <- list(
        Endpoints = list(
            Normalize_to_Coverage = list(
                type = "checkboxInput",
                value = TRUE
            ),
            Min_Normalization_Coverage = list(
                type = "numericInput",
                value = 100
            )
        )
    )
    getSvEndpoints <- function(conditions, groupLabels, targets){
        svs <- svCapture_regroupedSvs(
            input,
            assemblyOptions,
            assembly, groupedProjectSamples, groupingCols,
            conditions, groupLabels
        )
        svs$countGroups <- TRUE
        svs$dt <- svs$dt[, 
            .(
                POS = c(POS_1, POS_2) / 1e6,
                coverage = 1
            ), 
            by = .(TARGET_REGION, groupLabel)
        ]
        svs
    }
    normalizeSvEndpoints <- function(x, conditions, groupLabels){
        normalize <- endpointsPlot$plot$settings$get("Endpoints","Normalize_to_Coverage")
        if(!normalize) return(x)
        minCoverage <- endpointsPlot$plot$settings$get("Endpoints","Min_Normalization_Coverage")
        targets <- svCapture_assemblyTargets(assembly)
        coverage <- svCapture_getTargetCoverages(
            assembly, groupedProjectSamples, groupingCols, "all samples together", NULL, 
            conditions, groupLabels
        )
        for(trackLabel_ in targets$trackLabel){
            cov <- coverage$dt$dt[track == trackLabel_]
            cov[, y := y * 1e3]
            setkey(cov, x)
            x$dt[trackLabel == trackLabel_, coverage := sapply(x, function(ex) {
                cov[x >= ex][1, if(y >= minCoverage) 1 / y else NA_real_]
            })]
        }
        x
    }
    endpointsPlot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "SV Endpoint Coordinate (Mbp)",
        eventPlural = "ends",
        insideWidth = 1.25, 
        insideHeightPerBlock = 0.625,
        trackCols = "trackLabel",
        trackSameXLim = FALSE,
        trackSameYLim = TRUE,
        xlim = svCapture_xlimByTarget(assembly, "endpoints", endpointsPlot),
        extraSettings = endpointSettings,
        defaultBinSize = 10000 / 1e6,
        v = svCapture_vLinesByTarget(assembly),
        vShade = svCapture_vShadeByTarget(assembly),
        vColor = "grey40",
        x0Line = FALSE,
        y0Line = TRUE,
        aggFn = sum,
        aggCol = "coverage",
        # ylab = "SV Frequency",
        fontSize = 6.5,
        trackLabelPosition = "center",
        dataFn = function(conditions, groupLabels) {
            svCapture_parseDataByTarget(
                input,
                assemblyOptions,
                assembly, groupedProjectSamples, groupingCols,
                conditions, groupLabels,
                xCol = "POS", 
                dataFn = getSvEndpoints
            ) %>% normalizeSvEndpoints(conditions, groupLabels)
        }
    )
}
#----------------------------------------------------------------------
# target region coverage profiles
# ----------------------------------------------------------------------
svCapture_coverageServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    coverageSettings <- list(
        Coverage = list(
            Aggregate_Coverage = list(
                type = "selectInput",
                choices = c("by group","all samples together"),
                value = "all samples together"
            )
        )
    )
    coveragePlot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "Coordinate (Mbp)", 
        eventPlural = "samples",
        insideWidth = 1.25, 
        insideHeightPerBlock = 0.625,
        trackSameXLim = FALSE, 
        trackSameYLim = TRUE,
        xlim = svCapture_xlimByTarget(assembly, "coverage"),
        extraSettings = coverageSettings,
        defaultBinSize = 100 / 1e6,
        v = svCapture_vLinesByTarget(assembly),
        vShade = svCapture_vShadeByTarget(assembly),
        vColor = "grey40",
        x0Line = FALSE,
        y0Line = TRUE,
        pregrouped = TRUE,
        ylab = "Coverage (K)",
        fontSize = 6.5,
        trackLabelPosition = "center",
        dataFn = function(...) svCapture_getTargetCoverages(
            assembly, groupedProjectSamples, groupingCols, NULL, coveragePlot$plot, ...
        )
    )
}
#----------------------------------------------------------------------
# SV arc plots
# ----------------------------------------------------------------------
svx_getArcsPlotMetadata <- function(svs, settings){
    maxOffset <- trimws(settings$get("Arcs","Max_Offset"))
    maxOffset <- as.integer(if(maxOffset == "auto" || maxOffset == "") svs$maxOffset else maxOffset)
    interTargetColors <- CONSTANTS$plotlyColors[1:4]
    names(interTargetColors) <- c("LL","RR","LR","RL")
    list(
        maxOffset   = maxOffset,
        lwd         = settings$get("Arcs","Arc_Line_Width"),
        alpha       = settings$get("Arcs","Arc_Color_Alpha"),
        maxLines    = settings$get("Arcs","Max_Lines_Per_Type"),
        minAlphaLines = 50, # need at least this many lines to add transparency,
        interTargetColors = interTargetColors
    )
}
svx_arcs_addTopLegend <- function(plot, jxnTypes){
    altCodes <- svx_jxnTypes$altCode[svx_jxnTypes$altCode %in% jxnTypes]
    plot$addLegend(
        legend = svx_jxnType_altCodeToX(altCodes, "longName"),
        col =    svx_jxnType_altCodeToX(altCodes, "color"),
        x = "top",
        horiz = TRUE,
        x.intersp = 0,
        lwd = 1.5,
        lty = 1
    )
}
svx_getArcsTitle <- function(assembly, settings, family = "Plot"){
    title <- settings$get(family, "Title", NA)
    if(!isTruthy(title) && !is.null(assembly)) title <- assembly$env$DATA_NAME
    if(isTruthy(title)) title else ""
}
# ----------------------------------------------------------------------
svx_plotArcs_intraTarget <- function(plot, svs, input, assembly = NULL){
    md <- svx_getArcsPlotMetadata(svs, plot$settings)
    if(is.null(svs$intraTarget)) return()
    x <- seq(0, pi, length.out = 25)
    arc <- 0
    arcs <- do.call(rbind, lapply(svx_jxnType_longNameToX(input$Data_Types, "altCode"), function(jxnType){
        upsideDown <- jxnType != "L" # all but deletions plotted below
        svs_ <- svs$intraTarget[JXN_TYPE == jxnType]
        nSvs <- nrow(svs_)
        if(nSvs == 0) return(NULL)
        nLines <- min(md$maxLines, nSvs)
        alpha <- if(nLines <= md$minAlphaLines) 1 
                 else if(nLines == md$maxLines) md$alpha
                 else md$alpha + (1 - md$alpha) * (nLines - md$minAlphaLines) / (md$maxLines - md$minAlphaLines)
        color <- addAlphaToColor(svx_jxnType_altCodeToX(jxnType, "color"), alpha)
        svs_[sample(nSvs, nLines), {
            halfsize <- (OFFSET_2 - OFFSET_1) / 2
            center <- OFFSET_1 + halfsize
            y <- sin(x) * halfsize
            if(upsideDown) y <- -y
            arc <<- arc + 1
            .(
                x = cos(x) * halfsize + center, 
                y = y,
                col = color,
                arc = arc,
                jxnType = jxnType
            ) 
        }, by = .(PROJECT, SV_ID)]
    }))
    maxY <- max(abs(arcs$y))
    ylim <- c(-maxY * 1.05, maxY * 1.15)
    plot$initializeFrame(
        xlim = c(-md$maxOffset, md$maxOffset),
        ylim = ylim,
        xlab = "Distance from Gene Center (bp)",
        ylab = "",
        title = svx_getArcsTitle(assembly, plot$settings),
        xaxs = "i",
        yaxs = "i",
        yaxt = "n"
    )  
    for(arc_ in sample(arc)) {
        arc <- arcs[arc == arc_]
        lines(
            x = arc$x, 
            y = arc$y,
            col = arc$col[1],
            lwd = md$lwd
        )
    }
    abline(h = 0)
    svx_arcs_addTopLegend(plot, unique(arcs$jxnType))
}
svx_plotArcs_interTarget <- function(plot, svs, input, assembly = NULL){
    md <- svx_getArcsPlotMetadata(svs, plot$settings)
    if(is.null(svs$interTarget)) return()
    svs_ <- svs$interTarget
    nSvs <- nrow(svs_)
    if(nSvs == 0) return()
    nLines <- min(md$maxLines, nSvs)
    alpha <- if(nLines <= md$minAlphaLines) 1 
             else if(nLines == md$maxLines) md$alpha
             else md$alpha + (1 - md$alpha) * (nLines - md$minAlphaLines) / (md$maxLines - md$minAlphaLines)
    plot$initializeFrame(
        xlim = c(-md$maxOffset, md$maxOffset),
        ylim = c(-md$maxOffset, md$maxOffset),
        xlab = "Distance from Gene Center (bp)",
        ylab = "Distance from Gene Center (bp)",
        title = svx_getArcsTitle(assembly, plot$settings),
        xaxs = "i",
        yaxs = "i"
    ) 
    svs_ <- svs_[sample(nSvs, nLines)]
    colors <- svs_[, mapply(function(s1, s2){
        addAlphaToColor(md$interTargetColors[[paste0(s1, s2)]], alpha)
    }, SIDE_1, SIDE_2)]
    points(svs_$OFFSET_1, svs_$OFFSET_2, pch = 19, cex = 1, col = colors)
}
#----------------------------------------------------------------------
svCapture_arcsServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    settings <- c(assemblyPlotFrameSettings, list(
        Arcs = list(
            Arc_Plot_Type = list(
                type = "selectInput",
                choices = c(
                    "intraTarget",
                    "interTarget"
                ),
                value = "intraTarget"
            ),
            Max_Offset = list(
                type = "textInput",
                value = "auto"
            ), 
            Arc_Line_Width = list(
                type = "numericInput",
                value = 1
            ),
            Arc_Color_Alpha = list(
                type = "numericInput",
                value = 0.35
            ),
            Max_Lines_Per_Type = list(
                type = "numericInput",
                value = 500
            )
        )
    ))
    mars <- list(
        intraTarget = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$nullMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        interTarget = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$nullMar, 
            CONSTANTS$assemblyPlots$nullMar
        )
    )
    dims <- list(
        intraTarget = list(
            width  = 3, 
            height = 1.5
        ),
        interTarget = list(
            width  = 2, 
            height = 2
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
            targets <- svCapture_assemblyTargets(assembly)
            svs <- svCapture_matchingAssemblySvs(
                input,
                assemblyOptions, 
                assembly, groupedProjectSamples, groupingCols
            ) %>% svCapture_invertArcsByGene(assembly, assemblyOptions, targets)
            list(
                intraTarget = if(is.null(svs$intraTarget)) NULL else svs$intraTarget[groupLabel %in% groupLabels],
                interTarget = if(is.null(svs$interTarget)) NULL else svs$interTarget[groupLabel %in% groupLabels],
                maxOffset = svs$maxOffset
            )
        },  
        plotFrameFn = function(data) {
            plotAs <- assemblyPlot$plot$settings$get("Arcs","Arc_Plot_Type")
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
                plotAs <- assemblyPlot$plot$settings$get("Arcs","Arc_Plot_Type")
                assembly <- assembly()
                req(d, d$data, plotAs, assembly)
                startSpinner(session, message = paste("rendering", id))
                par(mar = plotFrameReactive()$mar)
                svs <- d$data
                get(paste("svx_plotArcs", plotAs, sep = "_"))(assemblyPlot$plot, svs, input, assembly)
                stopSpinner(session)
            }
        )
    )
}
#----------------------------------------------------------------------
# junction base usage profiles
# ----------------------------------------------------------------------
svCapture_junctionBasesSettings <- list(
    Junction_Bases = list(
        Max_Plotted_Distance = list(
            type = "numericInput",
            value = 25
        ),
        Distance_Grid_Spacing = list(
            type = "numericInput",
            value = 10
        ),
        Microhomology_Type = list(
            type = "selectInput",
            choices = c(
                "all",
                "microhomology",
                "insertion",
                "blunt",
                "no_insertions"
            ),
            value = "no_insertions"
        ),
        Transcription_Direction = list(
            type = "selectInput",
            choices = c(
                "all",
                "co_directional",
                "collisional"
            ),
            value = "all"
        ),
        Normalize_Base_Groups = list(
            type = "checkboxInput",
            value = FALSE
        )
    )
)
svCapture_junctionBasesSVs <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups,
    plot, conditions, groupLabels
){
    
    # filter SVs appropriate for profiling junction base usage, apply user microhomology type filter
    uHomType <- plot$settings$get("Junction_Bases", "Microhomology_Type")
    getMicrohomologyType <- function(sv){
        if(sv$MICROHOM_LEN > 0) "microhomology" 
        else if(sv$MICROHOM_LEN < 0) "insertion" 
        else "blunt"
    }
    svs <- svCapture_regroupedSvs(
        input,
        assemblyOptions,
        assembly, groupedProjectSamples, groupingCols,
        conditions, groupLabels
    )$dt[, 
        microhomologyType := getMicrohomologyType(.SD), by = .(PROJECT, SV_ID)
    ][switch(
        uHomType,
        all = TRUE,
        no_insertions = microhomologyType != "insertion",
        microhomologyType == uHomType
    )]  

    # determine SV gene strands
    targets <- copy(svCapture_assemblyTargets(assembly))
    setkey(targets, name) 
    svs[, 
        ":="(
            svKey = paste(project, SV_ID, sep = "::"),
            geneStrand = targets[TARGET_REGION, geneStrand]
        )
    ]
    setkey(svs, svKey)  
    svs
}
svCapture_filterFeatureMatrix <- function(featureMatrix, plot, svs){
    featureSvKeys <- substr(rownames(featureMatrix), 1, sapply(gregexpr("::", rownames(featureMatrix)), function(x) x[2] - 1)) # all SVs with a profile
    I <- featureSvKeys %in% svs$svKey

    # if requested, restrict the plot to one gene transcription orientation
    txnDirection <- plot$settings$get("Junction_Bases", "Transcription_Direction")
    if(txnDirection != "all"){
        featureMatrix <- featureMatrix[I, ]
        featureSvKeys <- featureSvKeys[I]
        x <- sapply(gregexpr("::", rownames(featureMatrix)), function(x) x[2] + 2)
        breakpointN <- as.integer(substr(rownames(featureMatrix), x, x)) 
        txnDirFn <- if(txnDirection == "co_directional"){ # fore each SV, keep the breakpoint where replication had the requested orientation relative to transcription
            function(geneStrand) if(geneStrand == "+") 1 else 2
        } else {
            function(geneStrand) if(geneStrand == "+") 2 else 1
        }   
        svs[, txnDirN := sapply(geneStrand, txnDirFn)]
        I <- svs[featureSvKeys, txnDirN == breakpointN]
    } else txnDirection <- ""
    list(
        matrix = featureMatrix[I, ],
        svKeys = featureSvKeys[I],
        txnDirection = txnDirection,
        nSvs = length(unique(featureSvKeys[I])),
        nBreakpoints = sum(I)
    )
}
svCapture_getBreakpointsTitle <- function(plot, d, env){
    title <- plot$settings$get("Plot", "Title", NA)
    nBreakpoints <- paste(commify(d$nBreakpoints), "breakpoints")
    if(!isTruthy(title)) title <- env$DATA_NAME
    title <- if(isTruthy(title)) paste0(title, " ", d$txnDirection, " (", nBreakpoints, ")")
                else nBreakpoints
}
svCapture_initJunctionPositionPlot <- function(plot, d, env, ylim, ylab){
    maxDist  <- plot$settings$get("Junction_Bases", "Max_Plotted_Distance")
    vSpacing <- plot$settings$get("Junction_Bases", "Distance_Grid_Spacing")
    plot$initializeFrame(
        xlim = c(-maxDist, maxDist),
        ylim = ylim,
        xlab = "Distance from Junction (bp)",
        ylab = ylab, 
        title = svCapture_getBreakpointsTitle(plot, d, env)
    )
    lastLine <- floor(maxDist / vSpacing) * vSpacing
    for(x in seq(-lastLine, lastLine, vSpacing)) lines(c(x, x), ylim, col = "grey80")
    lines(c(0, 0), ylim, col = "grey50")
    maxDist
}
svCapture_junctionBasesTopLegend <- function(plot, legend, col, lwd){
    plot$addLegend(
        legend = legend,
        col =    col,
        x = "top",
        horiz = TRUE,
        # x.intersp = 0,
        lwd = lwd,
        lty = 1
    )
}
# ----------------------------------------------------------------------
svCapture_junctionBasesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    getBaseUsage <- function(plot, conditions, groupLabels){ # dataFn
        startSpinner(session, message = "parsing base usage")
        svs <- svCapture_junctionBasesSVs(
            id, session, input, output, 
            isProcessingData, assemblyOptions,
            sourceId, assembly, groupedProjectSamples, groupingCols, groups,
            plot, conditions, groupLabels
        )

        # collect the base usage profile of the matching SVs
        d <- svCapture_filterFeatureMatrix(assembly()$baseUsageProfile, plot, svs) # [SV_side, position], value = numeric base, rownames = paste(project, svId, breakpointN)

        # aggregate all SV breakpoints along the span from retained to lost
        # bases = int [1:14298, 1:500] 1 2 2 3 4 3 2 3 4 1 ..., rowNames =  chr [1:14298] "HCT_ART558_052422::1496:1::1" "HCT_ART558_052422::1496:1::2" "HCT_ART558_052422::1503:1::1" "HCT_ART558_052422::1503:1::2" ...
        # counts = int [1:4, 1:500] 4296 2860 2870 4272 4339 2880 2901 4178 4350 2905 ... # [numeric base, position], value = number of SVs matching the row's base at position
        # fractions = num [1:4, 1:500] 0.3 0.2 0.201 0.299 0.303 ... # [numeric base, position], value = fraction of SVs matching the row's base at position
        # medians = num [1:4] 0.299 0.202 0.203 0.296 # [numeric base], value = median fraction by base
        counts <- apply(d$matrix, 2, function(svBasesAtPosition) sapply(1:4, function(base) sum(svBasesAtPosition == base, na.rm = TRUE))) 
        fractions <- t(t(counts) / colSums(counts, na.rm = TRUE))
        medians <- apply(fractions, 1, median, na.rm = TRUE)
        c(
            d[c("txnDirection","nSvs","nBreakpoints")], 
            list(
                baseUsageProfile = list(
                    counts = counts,
                    fractions = fractions,
                    medians = medians,
                    normalizedFractions = fractions / medians
                )
            )
        )

        # # summarize the base content in junction microhomologies
        # # these bases ARE already plotted as they are part of the retained sequence span
        # jxnBaseUsage <- assembly_$junctionBaseUsageProfile[, # c("PROJECT", "SV_ID", "microhomologyType", "A", "C", "G", "T")
        #     jxnSvKey := paste(PROJECT, SV_ID, sep = "::")
        # ][
        #     jxnSvKey %in% basesSvKeys[I] & 
        #     microhomologyType == "microhomology"
        # ][,
        #     geneStrand := svs[jxnSvKey, geneStrand]
        # ]
        # jxnBases <- jxnBaseUsage[,
        #     .(
        #         A  = ifelse(geneStrand == "+", A, T), # thus, microhomology bases are now all reported in transcription orientation
        #         C  = ifelse(geneStrand == "+", C, G), 
        #         G  = ifelse(geneStrand == "+", G, C), 
        #         T  = ifelse(geneStrand == "+", T, A),
        #         A_ = ifelse(geneStrand == "+", T, A), # while prime (_) are complemented
        #         C_ = ifelse(geneStrand == "+", G, C), 
        #         G_ = ifelse(geneStrand == "+", C, G), 
        #         T_ = ifelse(geneStrand == "+", A, T)
        #     )
        # ]
        # jxnBasesNames <- names(jxnBases)
        # jxnBases <- apply(jxnBases, 2, sum)
        # names(jxnBases) <- jxnBasesNames

        # # determine the position at the center of the microhomology base weighting
        # microhomologyPos <- svs[unique(jxnBaseUsage$jxnSvKey), .(mPos = 1:MICROHOM_LEN), by = .(svKey)][, mean(mPos)]
    }
    plotBaseUsage <- function(plot, d){ # plotFn
        d <- d$data
        assembly <- assembly()
        baseColors <- c(
            A = CONSTANTS$plotlyColors$red,
            C = CONSTANTS$plotlyColors$blue,
            G = CONSTANTS$plotlyColors$orange,
            T = CONSTANTS$plotlyColors$green
        )
        normalize <- plot$settings$get("Junction_Bases", "Normalize_Base_Groups")

        valueKey <- if(normalize) "normalizedFractions" else "fractions"
        values <- d$baseUsageProfile[[valueKey]]
        ylim <- range(values)
        ylim <- c(ylim[1], ylim[2] + diff(ylim) / 4)
        maxDist <- svCapture_initJunctionPositionPlot(
            plot, d, assembly$env, ylim, 
            if(normalize) "Normalized Base Fraction" else "Base Fraction"
        )

        lwd = 1
        if(!normalize) for(baseI in 1:4) plot$addLines(
            x = c(-maxDist, maxDist),
            y = rep(d$baseUsageProfile$medians[baseI], 2),
            col = baseColors[baseI],
            lwd = lwd
        ) else abline(h = 1, col = "grey50")

        x <- 1:ncol(values) - assembly$env$BASE_USAGE_SPAN
        lwd <- 1.5
        for(baseI in 1:4) plot$addLines(
            x = x,
            y = values[baseI, ],
            col = baseColors[baseI],
            lwd = lwd
        )
        baseOrder <- c("A", "T", "G", "C")        
        svCapture_junctionBasesTopLegend(plot, baseOrder, baseColors[baseOrder], lwd)

        # if(!normalize && d$txnDirection != ""){
        #     jxnBaseVals <- if(d$txnDirection == "co_directional") c("A","C","G","T") else c("A_","C_","G_","T_")
        #     jxnBaseVals <- d$jxnBases[jxnBaseVals]
        #     plot$addPoints(
        #         x = rep(-d$microhomologyPos, 4),
        #         y =  jxnBaseVals/ sum(jxnBaseVals),
        #         col = baseColors,
        #         pch = 1,
        #         cex = 1.5
        #     )
        # }


        # plot$addLegend(
        #     legend = baseOrder,
        #     col =    baseColors[baseOrder],
        #     x = "top",
        #     horiz = TRUE,
        #     x.intersp = 0,
        #     lwd = lwd,
        #     lty = 1
        # )
    }
    baseUsagePlot <- assemblyXYPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        extraSettings = svCapture_junctionBasesSettings,
        dims = list(width = 1.5, height = 1.5),
        mar = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        dataFn = getBaseUsage, 
        plotFn = plotBaseUsage
    )
}
#----------------------------------------------------------------------
# junction genome feature profiles
# ----------------------------------------------------------------------
svCapture_genomeFeaturesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    genomeFeatureSettings <- svCapture_junctionBasesSettings
    genomeFeatureSettings$Junction_Bases$Normalize_Base_Groups <- NULL
    genomeFeatureSettings <- c(genomeFeatureSettings, list(
        Genome_Features = list(
            Feature_Type = list(
                type = "selectInput",
                choices = c(
                    "flexibility",
                    "user_features"
                ),
                value = "flexibility"
            ),
            User_Feature = list(
                type = "textInput",
                value = ""
            )
        )
    ))
    getGenomeFeatures <- function(plot, conditions, groupLabels){ # dataFn
        startSpinner(session, message = "parsing genome features")
        svs <- svCapture_junctionBasesSVs(
            id, session, input, output, 
            isProcessingData, assemblyOptions,
            sourceId, assembly, groupedProjectSamples, groupingCols, groups,
            plot, conditions, groupLabels
        )

        # pick the dataset(s) to profile
        assembly_ <- assembly()
        featureType <- plot$settings$get("Genome_Features", "Feature_Type")
        featureMatrices <- if(featureType == "flexibility") list(
            flexibility = assembly_$flexibilityProfile
        ) else {
             # G-quadraplex, low-complexity 
            userFeature <- trimws(plot$settings$get("Genome_Features", "User_Feature"))
            if(!isTruthy(userFeature) || userFeature == "" || userFeature == "auto" || userFeature == "all"){
                userFeature <- names(assembly_$genomeFeatures)
            } else {
                req(userFeature %in% names(assembly_$genomeFeatures))
            }
            assembly_$genomeFeatures[userFeature]
        }

        # collect the feature profile(s) of the matching SVs
        d <- sapply(featureMatrices, svCapture_filterFeatureMatrix, plot, svs, simplify = FALSE, USE.NAMES = TRUE)
        sapply(d, function(dd){
            y <- apply(dd$matrix, 2, mean, na.rm = TRUE)
            c(
                dd[c("txnDirection","nSvs","nBreakpoints")], 
                list(
                    y = y,
                    minY = min(y),
                    maxY = max(y),
                    median = median(y)
                )
            )
        }, simplify = FALSE, USE.NAMES = TRUE)
    }
    plotGenomeFeatures <- function(plot, d){ # plotFn
        d <- d$data
        assembly <- assembly()
        maxDist  <- plot$settings$get("Junction_Bases", "Max_Plotted_Distance")
        vSpacing <- plot$settings$get("Junction_Bases", "Distance_Grid_Spacing")

        ylim <- range(sapply(d, function(dd) c(dd$minY, dd$maxY)))
        ylim <- c(ylim[1], ylim[2] + diff(ylim) / 4)
        featureType <- plot$settings$get("Genome_Features", "Feature_Type")
        maxDist <- svCapture_initJunctionPositionPlot(
            plot, d[[1]], assembly$env, ylim, 
            if(featureType == "flexibility") "Mean Flexibility Score" else "Fraction of SVs"
        )

        traceColors <- unlist(CONSTANTS$plotlyColors[c("blue","green","orange","red","black")])
        lwd = 1
        sapply(1:length(d), function(i){
            plot$addLines(
                x = c(-maxDist, maxDist),
                y = rep(d[[i]]$median, 2),
                col = traceColors[i],
                lwd = lwd
            )
        })
        lwd <- 1.5
        sapply(1:length(d), function(i){
            plot$addLines(
                x = 1:length(d[[i]]$y) - assembly$env$FEATURES_SPAN,
                y = d[[i]]$y,
                col = traceColors[i],
                lwd = lwd
            )
        })
        svCapture_junctionBasesTopLegend(plot, names(d), traceColors[1:length(d)], lwd)
    }
    genomeFeaturesPlot <- assemblyXYPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        extraSettings = genomeFeatureSettings,
        dims = list(width = 1.5, height = 1.5),
        mar = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        dataFn = getGenomeFeatures, 
        plotFn = plotGenomeFeatures
    )
}
#----------------------------------------------------------------------
# insertion template locations and properties
# ----------------------------------------------------------------------
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
    typePlotOrders <- list(
        foldback     = 4L, 
        palindrome   = 5L,
        crossJxn     = 3L,
        slippage     = 6L,
        strandSwitch = 2L, 
        other        = 1L  
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
        pileupSort       = svx_getInsertionsProperty(NULL, settings, NULL, "Pileup_Sort", default = "templateDistance"),
        flankUHom = sapply(1:maxInsertionSize, function(insertionSize){
            if(insertionSize >= minTemplateSize) 1L
            else as.integer(ceiling((minTemplateSize - insertionSize) / 2))
        }),
        nCombinations = sapply(1:maxInsertionSize, function(insertionSize){
            flankUHom <- if(insertionSize >= minTemplateSize) 1L
                        else as.integer(ceiling((minTemplateSize - insertionSize) / 2))
            4 ** (insertionSize + 2 * flankUHom)
        }),
        typeColors  = typeColors,
        typeLabels  = typeLabels,
        typeStrands = typeStrands,
        typePlotOrders = typePlotOrders
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
        xlab = "Insert Size (bp)",
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
        y = yield$trialSuccessProb * 100,
        pch = 19,
        col = "white"
    )
    plot$addPoints(
        x = -yield$MICROHOM_LEN,
        y = yield$trialSuccessProb * 100,
        pch = rev(as.character(md$flankUHom)),
        col = CONSTANTS$plotlyColors$black,
        cex = 0.85
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
        lines(c(x, x), c(0, 4.75), col = "grey80")
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
        col = sapply(svs$templateType, function(tt) md$typeColors[[tt]]),
        cex = 0.75
    )
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
            templateDistance,
            templateIsRc,
            series = paste(templateBreakpointN, templateType, sep = ":"),
            startPos = templateStartPos - md$searchSpace - if(templateBreakpointN == 1) 0L else 1L,
            endPos   = templateEndPos   - md$searchSpace - if(templateBreakpointN == 1) 0L else 1L,
            size     = templateEndPos - templateStartPos + 1,
            plotOrder = md$typePlotOrders[[templateType]]
        ), 
        by = .(PROJECT, SV_ID)
    ][, ":="(
        insStartPos = startPos + templateStartFm,
        insEndPos   = endPos   - templateEndFm
    )][
        abs(startPos) <= md$maxDist | 
        abs(endPos)   <= md$maxDist
    ]
    svs <- if(md$pileupSort == "size") svs[order(-size)] else svs[order(templateDistance, -size)]
    svs <- svs[,
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
    abline(v = seq(-1000, 1000, md$vSpacing), col = "grey80")
    abline(v = 0, col = CONSTANTS$plotlyColors$black)
    abline(h = c(1,3) * maxNSvs, col = CONSTANTS$plotlyColors$grey)
    abline(h = 2 * maxNSvs, col = CONSTANTS$plotlyColors$black)

    templateTypeCounts <- function(tt, bpn, isRc) {
        svs[
            templateType == tt & 
            templateBreakpointN == bpn &
            templateIsRc == isRc,
            if(.N == 0) NULL else paste(.N, md$typeLabels[[tt]], "   ")
        ]
    }
    if(md$maxDist < 50) for(bpn in 1:2) for(str in c(-1, 1)) text(
        if(bpn == 1) md$maxDist / 4 else -md$maxDist / 4 * 3, 
        (0.5 + (if(str == 1) 1 else 0) + (if(bpn == 1) 2 else 0)) * maxNSvs,
        paste(
            c(
                if(str == -1) character() else if(bpn == 1) "Gene Proximal Breakpoint" else "Gene Distal Breakpoint",
                sapply(templateTypes[md$typeStrands == str], templateTypeCounts, bpn, str == -1)
            ),
            collapse = "\n"
        ),
        adj = 0
    )
    # text( md$maxDist / 2, 3.6 * maxNSvs, "Gene Proximal\nBreakpoint")
    # text(-md$maxDist / 2, 0.6 * maxNSvs, "Gene Distal\nBreakpoint")
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
    }, keyby = .(plotOrder, series)]
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
    plot$initializeFrame(
        xlim = c(-md$maxDist, md$maxDist),
        ylim = c(0, d[, max(.SD, na.rm = TRUE) * 1.5, .SDcols = templateTypes]),
        xlab = "Distance from Junction (bp)",
        ylab = "# Templates",
        # yaxt = "n",
        bty = "n",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(nrow(svs)))
    )
    # abline(v = seq(vSpacing, maxDist - 1, vSpacing), col = CONSTANTS$plotlyColors$grey)
    abline(v = 0, col = CONSTANTS$plotlyColors$black)
    lwd <- 1.5 # plot$settings$Points_and_Lines()$Line_Width$value
    for(bpn in 1:2) for(type in templateTypes){
        dd <- d[templateBreakpointN == bpn]
        lines(dd$pos, dd[[type]], col = md$typeColors[[type]], lwd = lwd, lty = bpn)
    }
    svx_insertions_addTopLegend(plot, md)
}
svx_plotInsertions_microhomologyLengths <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    svs <- svs[templateType != "notFound"]
    xylim <- c(0.4, 8.6)
    plot$initializeFrame(
        xlim = xylim,
        ylim = xylim,
        xlab = "Priming Microhomology (bp)",
        ylab = "Resolving Microhomology (bp)",
        title = svx_getInsertionsTitle(assembly, plot$settings, nrow(svs)),
        xaxs = "i",
        yaxs = "i"
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
    d <- svs[, .(
        x = switch(
            templateType,
            foldback     = if(templateBreakpointN == 1) templateEndFm else templateStartFm,
            palindrome   = NA_integer_, # these microhomology lengths are untrustworthy
            crossJxn     = if(templateBreakpointN == 1) templateEndFm else templateStartFm,
            slippage     = if(templateBreakpointN == 1) templateStartFm else templateEndFm,
            strandSwitch = if(templateBreakpointN == 1) templateEndFm else templateStartFm,
            other        = NA_integer_
        ),
        y = switch(
            templateType,
            foldback     = if(templateBreakpointN == 2) templateEndFm else templateStartFm,
            palindrome   = NA_integer_, # these microhomology lengths are untrustworthy
            crossJxn     = if(templateBreakpointN == 2) templateEndFm else templateStartFm,
            slippage     = if(templateBreakpointN == 2) templateStartFm else templateEndFm,
            strandSwitch = if(templateBreakpointN == 2) templateEndFm else templateStartFm,
            other        = NA_integer_
        ),
        templateType
    ), by = .(templateType, templateBreakpointN)]
    mdiLevelPlot(
        d,
        xlim = xylim,
        xinc = 1,
        ylim = xylim,
        yinc = 1,
        z.fn = function(d) sum(!is.na(d)),   # function applied to z.columnumn, per grid spot, to generate the output color
        z.column = "x", # the column in dt passed to z.fn, per grid spot
        settings = levelPlotSettings, # a settings object from the enclosing staticPlotBox, or any list compatible with mdiLevelPlotSettings
        legendTitle = "# of SVs",
        h = c(mode(d$y), quantile(d$y, c(0.5, 0.95), na.rm = TRUE)),
        v = c(mode(d$x), quantile(d$x, c(0.5, 0.95), na.rm = TRUE))
    )
}
svx_plotInsertions_templateSizes <- function(plot, svs, assembly = NULL){
    md <- svx_getInsertionPlotMetadata(assembly, plot$settings)
    svs <- svs[templateType != "notFound"]
    xlim <- c(1 - 0.6, md$maxInsertionSize + 4.6) #m d$minInsertionSize - 
    ylim <- c(1 - 0.6, md$maxInsertionSize + 4.6)  # md$minTemplateSize - 
    plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = "Insertion Size (bp)",
        ylab = "Total Template Size (bp)",
        title = svx_getInsertionsTitle(assembly, plot$settings, nrow(svs)),
        xaxs = "i",
        yaxs = "i"
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
    d <- svs[, .(
        x = -MICROHOM_LEN,
        y = -MICROHOM_LEN + templateStartFm + templateEndFm,
        templateType
    )]
    mdiLevelPlot(
        d,
        xlim = xlim,
        xinc = 1,
        ylim = ylim,
        yinc = 1,
        z.fn = function(d) sum(!is.na(d)),   # function applied to z.columnumn, per grid spot, to generate the output color
        z.column = "x", # the column in dt passed to z.fn, per grid spot
        settings = levelPlotSettings, # a settings object from the enclosing staticPlotBox, or any list compatible with mdiLevelPlotSettings
        legendTitle = "# of SVs",
        h = c(mode(d$y), quantile(d$y, c(0.5, 0.95), na.rm = TRUE)),
        v = c(mode(d$x), quantile(d$x, c(0.5, 0.95), na.rm = TRUE))
    )
}
svCapture_insertionTemplatesServer <- function(
    id, session, input, output, 
    isProcessingData, assemblyOptions,
    sourceId, assembly, groupedProjectSamples, groupingCols, groups
){
    settings <- c(assemblyPlotFrameSettings, list(
        Insertions = list(
            Plot_Insertions_As = list(
                type = "selectInput",
                choices = c(
                    "nTemplateInstances",
                    "yield",
                    "map",
                    "pileup",
                    "histogram",
                    "microhomologyLengths",
                    "templateSizes"
                ),
                value = "pileup"
            ),
            P_Value_Threshold = list(
                type = "numericInput",
                value = 0.05
            ),
            Max_Plotted_Distance = list(
                type = "numericInput",
                value = 25
            ),
            Distance_Grid_Spacing = list(
                type = "numericInput",
                value = 5
            ),
            Pileup_Sort = list(
                type = "selectInput",
                choices = c(
                    "size",
                    "templateDistance"
                ),
                value = "size"
            ),
            Spacer = list(
                type = "spacer"
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
        nTemplateInstances = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
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
        pileup = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$titleMar, 
            CONSTANTS$assemblyPlots$nullMar
        ),
        histogram = c(
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
        templateSizes = c(
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar, 
            CONSTANTS$assemblyPlots$stdAxisMar
        )
    )
    dims <- list(
        nTemplateInstances = list(
            width  = 1.5, 
            height = 1.25
        ),
        yield = list(
            width  = 1.5, 
            height = 1.25
        ),
        map = list(
            width  = 4, 
            height = 2
        ),
        pileup = list(
            width  = 4, 
            height = 3
        ),        
        histogram = list(
            width  = 4, 
            height = 1.5
        ),
        microhomologyLengths = list(
            width  = 1, 
            height = 1
        ),
        templateSizes = list(
            width  = 1, 
            height = 1
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
            targets <- svCapture_assemblyTargets(assembly)
            svs <- svCapture_matchingAssemblySvs(
                input,
                assemblyOptions, 
                assembly, groupedProjectSamples, groupingCols
            ) %>% svCapture_invertInsertionsByGene(assembly, assemblyOptions, targets)
            svs[
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