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
    xCol, expandFn = function(dt) dt
){
    targets <- svCapture_assemblyTargets(assembly)
    keyedTargets <- copy(targets)
    setkey(keyedTargets, name)
    svs <- svCapture_matchingAssemblySvs(
        input,
        assemblyOptions, 
        assembly, groupedProjectSamples, groupingCols
    ) %>% expandFn()
    dt <- do.call(rbind, lapply(targets$name, function(target){
        x <- regroupToUserConditions(
            svs[TARGET_REGION == target], 
            groupingCols, conditions, groupLabels
        )
        x$dt[, ":="(
            x = x$dt[[xCol]],
            trackLabel = keyedTargets[target, trackLabel]
        )]
    }))
    list(
        titleSuffix = NULL,
        groupLabels = unique(dt$groupLabel),
        groupCounts = dt[, .N, by = .(groupLabel)],
        trackLabels = targets$trackLabel,
        dt = dt
    )
}
svCapture_xlimByTarget <- function(assembly){
    reactive({
        targets <- svCapture_assemblyTargets(assembly)
        x <- lapply(1:nrow(targets), function(i) targets[i, c(paddedStart + 1, paddedEnd) / 1e6])
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
    plot <- assemblyDensityPlotServer(
        id, session, input, output, 
        isProcessingData, assemblyOptions,
        sourceId, assembly, groupedProjectSamples, groupingCols, groups,
        xlab = "SV Endpoint Coordinate (Mbp)",
        eventPlural = "SVs",        
        insideWidth = 2.5, 
        insideHeightPerBlock = 0.75,
        trackCols = "trackLabel",         
        trackSameXLim = FALSE,
        trackSameYLim = TRUE,
        xlim = svCapture_xlimByTarget(assembly),        
        defaultBinSize = 10000 / 1e6,
        v = svCapture_vLinesByTarget(assembly),
        vShade = svCapture_vShadeByTarget(assembly),
        vColor = "grey40",
        x0Line = FALSE,
        y0Line = TRUE,
        dataFn = function(conditions, groupLabels) {
            svCapture_parseDataByTarget(
                input,
                assemblyOptions,
                assembly, groupedProjectSamples, groupingCols,
                conditions, groupLabels,
                xCol = "POS", 
                expandFn = function(dt) dt[
                    groupLabel %in% groupLabels, 
                    .(POS = c(POS_1, POS_2) / 1e6) , 
                    by = .(TARGET_REGION, groupLabel)
                ]
            )
        }
    )
}
#----------------------------------------------------------------------
