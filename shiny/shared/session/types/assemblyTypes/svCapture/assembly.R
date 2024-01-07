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
    x <- svCapture_matchingAssemblySvs(
        input,
        assemblyOptions,
        assembly, groupedProjectSamples, groupingCols
    )[
        groupLabel %in% groupLabels
    ] %>% regroupToUserConditions(
        groupingCols, conditions, groupLabels
    )
    x$dt <- expandFn(x$dt)
    x$dt[, ":="(
        x = x$dt[[xCol]],
        trackLabel = unlist(keyedTargets[TARGET_REGION, trackLabel])
    )]
    x$trackLabels <- targets$trackLabel
    x$groupCounts <- x$dt[, .N, by = .(groupLabel)]
    x$trackCounts <- x$dt[, .N, by = .(trackLabel)]
    x 
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
        eventPlural = "ends",
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
                expandFn = function(dt) dt[, 
                    .(POS = c(POS_1, POS_2) / 1e6), 
                    by = .(TARGET_REGION, groupLabel)
                ]
            )
        }
    )
}
#----------------------------------------------------------------------
svx_getInsertionsProperty <- function(assembly, settings, assemblyKey, settingsKey, default = NA){
    if(is.null(assembly) || is.null(assembly$env[[assemblyKey]])){
        settingsFamily <- if("Junction_Properties" %in% names(settings)) "Junction_Properties" else "Insertions"
        if(is.null(settings[[settingsFamily]])) return(default)
        settings[[settingsFamily]]()[[settingsKey]]$value
    } else {
        assembly$env[[assemblyKey]]
    }
}
svx_getInsertionsTitle <- function(assembly, settings, nSvs, family = "Plot"){
    title <- settings$get("Plot", "Title", NA)
    nSvs <- paste(nSvs, "SVs")
    if(!isTruthy(title) && !is.null(assembly)) title <- assembly$env$DATA_NAME
    if(isTruthy(title)) paste0(title, " (", nSvs, ")")
    else nSvs
}
svx_plotInsertions_yield <- function(plot, svs, assembly = NULL){
    flankUHom   <- svx_getInsertionsProperty(assembly, plot$settings, "FLANKING_MICROHOMOLOGY", "Flanking_Microhomology", default = 2)
    searchSpace <- svx_getInsertionsProperty(assembly, plot$settings, "INSERTION_SEARCH_SPACE", "Insertion_Search_Space", default = 500)
    searchSpace <- searchSpace * 2 * 2 * 2 # two breakpoint junctions searched on both sides of the breakpoint on two strands 
    yield <- svs[, .(
        nCombinations = 4 ** (-MICROHOM_LEN + flankUHom * 2), # number of possible unique template sequences for an insertion size 
        nSvs = .N,                                # number of trials ...
        nFound = sum(templateType != "not_found") # number of successes ...
    ), keyby = MICROHOM_LEN][,                    # ... by insert size
        mu := searchSpace * (1 / nCombinations)   # Poisson mean, i.e., expected target density ...
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
    pThreshold <- 0.05
    plot$addPoints(
        x = -yield$MICROHOM_LEN,            
        y = yield$successRate * 100,
        col = ifelse(yield$pValue < 0.05, 
                        CONSTANTS$plotlyColors$red, 
                        CONSTANTS$plotlyColors$black)
    )
    plot$addLegend(
        x = "topright",
        legend = c(
            paste("p", "<",  pThreshold),
            paste("p", ">=", pThreshold),
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
    svs <- svs[templateType != "not_found"]
    svs[, templateType2 := templateType]
    svs[templateBreakpointSide == "lost", templateType2 := "other"]

    maxDist   <- svx_getInsertionsProperty(NULL, plot$settings, NULL, "Max_Plotted_Distance",  default = 2)
    vSpacing  <- svx_getInsertionsProperty(NULL, plot$settings, NULL, "Distance_Grid_Spacing", default = 5)
    flankUHom <- svx_getInsertionsProperty(assembly, plot$settings, "FLANKING_MICROHOMOLOGY", "Flanking_Microhomology", default = 2)
    foldback <- CONSTANTS$plotlyColors$blue
    crossJxn <- CONSTANTS$plotlyColors$green
    other    <- CONSTANTS$plotlyColors$black
    cols <- list(foldback = foldback, crossJxn = crossJxn, other = other)

    svs[, ":="(
        x = switch(
            templateBreakpointN,
            switch(templateBreakpointSide, lost =  templateDistance, retained = -templateDistance),
            switch(templateBreakpointSide, lost = -templateDistance, retained =  templateDistance)
        ),
        y = switch(
            templateBreakpointN,
            switch(templateType, crossJxn = 4, foldback = 3),
            switch(templateType, crossJxn = 2, foldback = 1)
        )
    ), by = c("PROJECT", "SV_ID")]

    plot$initializeFrame(
        xlim = c(-maxDist, maxDist),          
        ylim = c(0.5, 6),
        xlab = "Insertion Template Distance from Junction (bp)",            
        ylab = "",
        yaxt = "n",
        bty = "n",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(nrow(svs)))
    )
    lastLine <- floor(maxDist / vSpacing) * vSpacing
    for(x in seq(-lastLine, lastLine, vSpacing)){
        lines(c(x, x), c(0, 4.75), col = CONSTANTS$plotlyColors$grey)
    }
    boldLwd <- 2.5
    for(y in 3:4){ # left side of junction, drawn on top two lines
        lines(c(-lastLine, 0), c(y, y), col = CONSTANTS$plotlyColors$black, lwd = boldLwd)  
        lines(c(0, lastLine),  c(y, y), col = CONSTANTS$plotlyColors$grey)   
    }
    for(y in 1:2){ # right side of junction, drawn on bottom two lines
        lines(c(-lastLine, 0), c(y, y), col = CONSTANTS$plotlyColors$grey)  
        lines(c(0, lastLine),  c(y, y), col = CONSTANTS$plotlyColors$black, lwd = boldLwd)   
    }        
    lines(c(0, 0), c(1, 4), col = CONSTANTS$plotlyColors$black, lwd = boldLwd)
    text(-maxDist, 5.25, "Left Side of SV",  pos = 4)
    text( maxDist, 5.25, "Right Side of SV", pos = 2)
    plot$addPoints(
        x = jitter(svs$x, amount = 0.25), # jitter to make all points visible      
        y = jitter(svs$y, amount = 0.25),
        col = unlist(cols[svs$templateType2])
    )
    plot$addLegend(
        legend = c("cross-junction", "foldback", "other"),
        col =    c(crossJxn,         foldback,    other),
        x = "top",
        horiz = TRUE,
        x.intersp = 0,
        pch = 19,
        pt.cex = 1 
    )
}
svx_plotInsertions_histogram <- function(plot, svs, assembly = NULL){
    svs <- svs[templateType != "not_found"]
    svs[templateBreakpointSide == "lost", templateType := "other"]

    maxDist   <- svx_getInsertionsProperty(NULL, plot$settings, NULL, "Max_Plotted_Distance",  default = 2)
    vSpacing  <- svx_getInsertionsProperty(NULL, plot$settings, NULL, "Distance_Grid_Spacing", default = 5)
    flankUHom <- svx_getInsertionsProperty(assembly, plot$settings, "FLANKING_MICROHOMOLOGY", "Flanking_Microhomology", default = 2)
    foldback <- CONSTANTS$plotlyColors$blue
    crossJxn <- CONSTANTS$plotlyColors$green
    other    <- CONSTANTS$plotlyColors$black
    cols <- list(foldback = foldback, crossJxn = crossJxn, other = other)

    d <- dcast(svs, templateDistance ~ templateType, fun.aggregate = length, fill = 0)
    d <- merge(data.table(templateDistance = 1:maxDist), d, by = "templateDistance", all.x = TRUE)
    d[is.na(d)] <- 0

    plotCols <- c("crossJxn", "foldback", "other")
    plot$initializeFrame(
        xlim = c(0, maxDist),          
        ylim = c(0, d[, max(.SD, na.rm = TRUE) * 1.05, .SDcols = plotCols]),
        xlab = "Distance from Junction (bp)",            
        ylab = "# of Templates",
        title = svx_getInsertionsTitle(assembly, plot$settings, sum(nrow(svs)))
    )
    abline(v = seq(vSpacing, maxDist - 1, vSpacing), col = CONSTANTS$plotlyColors$grey)
    lwd <- 1.5 # plot$settings$Points_and_Lines()$Line_Width$value
    for(type in plotCols){
        lines(d$templateDistance, d[[type]], col = cols[[type]], lwd = lwd)
    }
    plot$addLegend(
        legend = c("cross-junction", "foldback", "other"),
        col =    unlist(cols[plotCols]),
        x = "topright",
        lty = 1,
        lwd = lwd
    )
}
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
                    "histogram"
                ),
                value = "map"
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
            width  = 1.75, 
            height = 1.25
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
                !is.na(templateType) # "not_found", "foldback", or "crossJxn"
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
                get(paste("svx_plotInsertions", plotAs, sep = "_"))(assemblyPlot$plot, d$data, assembly)
                stopSpinner(session)
            }
        )
    )
}
