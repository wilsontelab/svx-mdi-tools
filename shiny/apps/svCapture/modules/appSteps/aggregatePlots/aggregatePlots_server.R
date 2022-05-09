#----------------------------------------------------------------------
# reactive components to plot summary results over all selected samples
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
aggregatePlotsServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'aggregatePlots' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    templates = list(
        file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), 
        file.path(app$sources$suiteGlobalDir, "settings", "svCapture_filters.yml"),
        id
    ),
    fade = FALSE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues()

#----------------------------------------------------------------------
# generate the list of all filtered SVs from all selected samples
#----------------------------------------------------------------------
targetClasses <- reactive({ SVX$targetClasses[[settings$SV_Filters()$Target_Class$value]] })
filteredSvs <- reactive({ 
    getFilteredSvs(settings, sampleSelector, isCapture = TRUE) 
})

#----------------------------------------------------------------------
# group and aggregate sv rate data
#----------------------------------------------------------------------
svRates <- reactive({

    # get samples and their categorical assignments
    samples <- sampleSelector$selectedSamples() 
    req(samples)       
    assignments <- sampleSelector$selectedAssignments()
    req(assignments)
    req(nrow(assignments) > 0)   

    # order by group/category
    assignments[, originalOrder := 1:.N]
    transpose <- settings$Data()$Transpose$value
    assignments <- if(transpose){
        assignments[order(Category2, Category1, originalOrder)]
    } else {
        assignments[order(Category1, Category2, originalOrder)]
    }

    # count the number of SVs assigned to each sample
    # depending on filters, an SV could be assigned to multiple samples
    ps <- filteredSvs()[, PROJECT_SAMPLES]
    assignments[, nSvs := sum(sapply(ps, function(x) uniqueId %in% x)), by = uniqueId]  

    # normalize SV counts to rates based on on-target coverage depth
    assignments[, onTargetCoverage := {
        manifest <- loadPersistentFile(sourceId = Source_ID, contentFileType = "manifestFile", sep = ",") 
        sid <- Sample_ID
        persistentCache[[manifest]]$data[Sample_ID == sid, onTargetCoverage]
    }, by = c("Source_ID", "Sample_ID")]
    assignments[, svRate := nSvs / onTargetCoverage]

    # assign replicate numbers
    sourceIds <- assignments[, unique(Source_ID)]
    replicates <- seq_along(sourceIds) 
    names(replicates) <- sourceIds   
    assignments[, replicate := replicates[Source_ID]]                    

    # simplify groups from the data to be plotted
    isCategory1 <- assignments[, length(unique(Category1)) > 1]
    isCategory2 <- assignments[, length(unique(Category2)) > 1]
    assignments[, group := if(isCategory1 && isCategory2){
        if(transpose){
            paste(Category2Name, Category1Name, sep = "\n")
        } else {
            paste(Category1Name, Category2Name, sep = "\n")
        }
    } else if(isCategory1){
        Category1Name
    } else if(isCategory2){
        Category2Name
    } else uniqueId]

    assignments
})

#----------------------------------------------------------------------
# SV rates plot
#----------------------------------------------------------------------
svRatesPlot <- staticPlotBoxServer(
    'svRates', 
    legend = TRUE,
    immediate = TRUE,
    points = TRUE,
    margins = TRUE,
    title = TRUE,
    create = function(...){
        svRates <- copy(svRates())[order(rev(1:.N))]
        req(svRates)
        svRates[, groupFactor := factor(group, unique(group))]
        xi <- as.integer(svRates$groupFactor)
        svRatesPlot$initializeFrame(
            xlim = c(0, max(svRates$svRate) * 1.1),            
            ylim = c(0, max(as.integer(svRates$groupFactor))) + 0.5,
            xlab = "SV Rate",            
            ylab = "",
            yaxt = "n"
        )
        abline(v = seq(0, 1, 0.1), col = "grey80")
        abline(h = c(0, xi) + 0.5, col = "grey80")
        svRatesPlot$addPoints(
            x = svRates$svRate,            
            y = jitter(xi, amount = 0.1),
            col = unlist(CONSTANTS$plotlyColors[as.integer(factor(svRates$replicate))])
        )
        axis(
            2, 
            at = xi, 
            labels = svRates$group,
            las = 1
        )
        # barplot(as.matrix(x), beside = TRUE, horiz = TRUE)
    }
)

#----------------------------------------------------------------------
# SV Size and Microhomology Length Distribution plots
#----------------------------------------------------------------------
svsByGroup <- reactive({
    svRates <- copy(svRates())
    req(svRates)
    setkey(svRates, uniqueId)  
    ps <- svRates$uniqueId
    filteredSvs()[, .(
        sequenced = JXN_BASES != "*",
        MICROHOM_LEN, 
        SV_SIZE,
        group = switch(
            settings$get("Data", "Distribute"),
            byGroup = {
                x <- unlist(PROJECT_SAMPLES)
                unique(svRates[x[x %in% ps], group])
            },
            bySample = {
                x <- unlist(PROJECT_SAMPLES)
                x[x %in% ps]
            },
            allSamplesTogether = "allSamplesTogether"
        )
    ), by = SV_ID]
})
plotFrequencyDistribution <- function(
    plot,
    column, 
    xlim,     
    xlab, 
    modifyData = function(...) NULL,  
    modifyPlotBase = function(...) NULL
){
    svsByGroup <- copy(svsByGroup()) # since we might modify it...
    req(svsByGroup)
    isCum <- settings$get("Data", "Cumulative")    
    modifyData(svsByGroup, isCum)
    filter <- if(column == "MICROHOM_LEN") svsByGroup$sequenced else TRUE

    d <- svsByGroup[filter, .N, keyby = c("group", column)]

    d[, Frequency := N / sum(N), keyby = group]
    ymax <- if(isCum) 1 else max(d$Frequency) * 1.05
    plot$initializeFrame(
        xlim = xlim,
        ylim = c(0, ymax),
        xlab = xlab,
        ylab = paste0(if(isCum) "Cumulative " else "", "Frequency")
    )
    modifyPlotBase(d)
    if(isCum) abline(h = 0.5, col = "grey60")
    groups <- unique(d$group)
    for(i in seq_along(groups)){
        dd <- d[group == groups[i]]
        plot$addLines(
            x = dd[[column]], 
            y = if(isCum) cumsum(dd$Frequency) else dd$Frequency, 
            col = CONSTANTS$plotlyColors[[i]]
        )
    }  
    if(settings$get("Data", "Distribute") != "allSamplesTogether"){
        plot$addLegend(
            legend = groups,
            col = unlist(CONSTANTS$plotlyColors[seq_along(groups)])
        )
    }
}
sizeDistribution <- staticPlotBoxServer(
    'sizeDistribution',
    legend    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){
        plotFrequencyDistribution(
            sizeDistribution,
            column  = "SV_SIZE", 
            xlim    = log10(c(1000, 1e6)),            
            xlab    = "log10 SV Size (bp)",
            modifyData = function(svsByGroup, isCum) {
                svsByGroup[, SV_SIZE := log10(SV_SIZE)]
                if(!isCum) svsByGroup[, SV_SIZE := round(SV_SIZE * 2, 0) / 2]
            },
            modifyPlotBase = function(...){
                abline(v = 0:10, col = "grey60")
            }
        )
    }
)
microhomologyDistribution <- staticPlotBoxServer(
    'microhomologyDistribution',
    legend    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){
        plotFrequencyDistribution(
            microhomologyDistribution,            
            column  = "MICROHOM_LEN", 
            xlim    = c(-20, 20),            
            xlab    = "Microhomology Length (bp)",
            modifyPlotBase = function(...){
                abline(v = seq(-100, 100, 10), col = "grey60")
                abline(v = 0)
            }
        )
    }
)

# ----------------------------------------------------------------------
# table of the plotted data
# ----------------------------------------------------------------------
ratesTable <- bufferedTableServer(
    id = 'ratesTable',
    parentId = id,
    parentInput = input,
    selection = 'single',
    tableData = reactive({
        svRates <- svRates()
        req(svRates)
        svRates[, .(
            Group       = Category1Name,
            Category    = Category2Name,
            Project     = Project,
            Sample      = Sample_ID,
            nSvs        = nSvs, 
            onTargetCoverage = onTargetCoverage,
            svRate      = svRate
        )]
    })
)

# ----------------------------------------------------------------------
# figures related to t- SVs
# ----------------------------------------------------------------------
distancePlot <- staticPlotBoxServer(
    'distancePlot',
    legend    = TRUE,
    points    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){

        # load targets bed, assumed to be the same for all sample sources
        assignments <- sampleSelector$selectedAssignments()
        req(assignments)
        req(nrow(assignments) > 0)
        sourceId <- assignments[, Source_ID[1]]
        targetsBed <- loadPersistentFile(sourceId = sourceId, contentFileType = "targetsBed") 
        targets <- persistentCache[[targetsBed]]$data
        targets[, center := start + size / 2]
        targetRegionNames <- c(targets$regionName, paste0("*,", targets$regionName))
        targetRegionCenters <- as.list(c(targets$center, targets$center))
        names(targetRegionCenters) <- targetRegionNames

        # TODO: expose as settings
        binSize <- 50000
        maxBins <- 100
        maxDistance <- binSize * maxBins

        # get the summed onTargetCoverage, i.e., TT Proper, for the libraries being analyzed, for normalization
        libraryMetrics <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
            manifestFile <- getSourceFilePath(sourceId, "manifestFile")
            x <- fread(manifestFile)
            x[, ':='(uniqueId = paste(Project, Sample_ID, sep = ":"))]
            x
        }))
        onTargetCoverage <- libraryMetrics[uniqueId %in% assignments[, uniqueId], sum(onTargetCoverage)]

        # load the relevant SVs
        svs <- getFilteredSvs(
            settings, 
            sampleSelector, 
            isCapture = TRUE, 
            targetClasses = SVX$targetClasses$distancePlot, # TT, TA, t- types
            noSize = TRUE
        ) 

        # calculate distance of the farthest SV junction point from the relevant target region center
        svs <- svs[
            CHROM_1 == CHROM_2 & # only del/dup/inv are informative here
            TARGET_REGION %in% targetRegionNames, 
        .(
            distanceBin = {
                distances <- c(POS_1, POS_2) - targetRegionCenters[[TARGET_REGION]]
                as.integer(distances[which.max(abs(distances))] / binSize) * binSize
            }
        ), by = SV_ID]

        # calculate the binned Tx SV count, normalized to TT Proper
        x <- svs[
            abs(distanceBin) <= maxDistance, 
            .(normalizedSvCount = .N / onTargetCoverage), 
            by = distanceBin
        ]
        # zeros <- data.table(distanceBin = seq(-maxDistance, maxDistance, binSize))
        x <- merge(
            x, 
            data.table(distanceBin = seq(-maxDistance, maxDistance, binSize)
        ), all.y = TRUE)[order(distanceBin)]
        x[is.na(normalizedSvCount), normalizedSvCount := 0]

        # construct the plot
        distancePlot$initializeFrame(
            xlim = c(-1, 1) * maxDistance, 
            ylim = c(0, max(x$normalizedSvCount)),
            xlab = "Distance of Farthest End from Target (bp)",
            ylab = "Normalized Bin SV Count"
        )
        abline(v = 0, col = CONSTANTS$plotlyColors$black)
        isZero <- x$normalizedSvCount == 0        
        distancePlot$addLines(
            x = x$distanceBin,
            y = x$normalizedSvCount,
            col = CONSTANTS$plotlyColors$grey
        )
        # distancePlot$addPoints(
        #     x = x$distanceBin[isZero],
        #     y = x$normalizedSvCount[isZero],
        #     col = CONSTANTS$plotlyColors$grey
        # )
        distancePlot$addPoints(
            x = x$distanceBin[!isZero],
            y = x$normalizedSvCount[!isZero],
            col = CONSTANTS$plotlyColors$blue
        )  
    }
)
ligationArtifactPlot <- staticPlotBoxServer(
    'ligationArtifactPlot',
    legend    = TRUE,
    points    = TRUE,
    lines     = TRUE,
    title     = TRUE,
    margins   = TRUE,
    immediate = TRUE,
    create = function(...){

        # load the relevant SVs
        svs <- getFilteredSvs(
            settings, 
            sampleSelector, 
            isCapture = TRUE, 
            targetClasses = "t-",
            noSize = TRUE
        ) 

        x <- svs[
            JXN_BASES != "*",
            .(freq = .N / nrow(svs)),
            by = MICROHOM_LEN         
        ]

        maxMH <- 20

        x <- merge(
            x, 
            data.table(MICROHOM_LEN = seq(-maxMH, maxMH, 1)
        ), all.y = TRUE)[order(MICROHOM_LEN)]
        x[is.na(freq), freq := 0]

        # construct the plot
        distancePlot$initializeFrame(
            xlim = c(-1, 1) * maxMH, 
            ylim = c(0, max(x$freq)) * 1.1,
            xlab = "Microhomology Length (bp)",
            ylab = "Frequency"
        )
        abline(v = 0, col = CONSTANTS$plotlyColors$black)
        # isZero <- x$normalizedSvCount == 0        
        distancePlot$addLines(
            x = x$MICROHOM_LEN,
            y = x$freq,
            col = CONSTANTS$plotlyColors$blue
        )
        # distancePlot$addPoints(
        #     x = x$distanceBin[isZero],
        #     y = x$normalizedSvCount[isZero],
        #     col = CONSTANTS$plotlyColors$grey
        # )
        # distancePlot$addPoints(
        #     x = x$distanceBin[!isZero],
        #     y = x$normalizedSvCount[!isZero],
        #     col = CONSTANTS$plotlyColors$blue
        # )  
    }
)

# ----------------------------------------------------------------------
# define bookmarking actions
# ----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    sampleSet <- bm$input[['sampleSelector-sampleSet']]
    sampleSelector$setSampleSet(sampleSet) 
    if(!is.null(bm$outcomes)) {
        outcomes <<- listToReactiveValues(bm$outcomes)
        sampleSelector$setSelectedSamples(sampleSet, bm$outcomes$samples)
        svRatesPlot$settings$replace(bm$outcomes$svRatesPlotSettings)
        sizeDistribution$settings$replace(bm$outcomes$sizeDistributionSettings)
        microhomologyDistribution$settings$replace(bm$outcomes$microhomologyDistributionSettings)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    samples  = sampleSelector$selectedSamples,
    outcomes = reactive({ list(
        samples = sampleSelector$selectedSamples(),
        svRatesPlotSettings = svRatesPlot$settings$all_(),
        sizeDistributionSettings = sizeDistribution$settings$all_(),
        microhomologyDistributionSettings = microhomologyDistribution$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
