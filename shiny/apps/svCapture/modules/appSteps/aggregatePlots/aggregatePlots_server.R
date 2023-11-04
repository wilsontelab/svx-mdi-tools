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
    fade = FALSE,
    presets = settingsPresets
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
    assignments <- copy(sampleSelector$selectedAssignments())
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
    uniqueIds <- svRates$uniqueId # all sv ids in plotted samples 
    groupType <- settings$get("Data", "Distribute")

    # expand SVs across samples or groups
    # a SV will have multiple rows in the output if it is assigned to multiple samples or groups
    expandedSvs <- filteredSvs()[, .( 
        sequenced = JXN_BASES != "*",
        MICROHOM_LEN, 
        SV_SIZE,
        group = switch(
            groupType,
            byGroup = {
                x <- unlist(PROJECT_SAMPLES)
                unique(svRates[x[x %in% uniqueIds], group])
            },
            bySample = {
                x <- unlist(PROJECT_SAMPLES)
                x[x %in% uniqueIds]
            },
            allSamplesTogether = "allSamplesTogether"
        )
    ), by = c("PROJECT", "SV_ID")]

    # add columns for sorting and naming in legends
    if(groupType == "byGroup"){
        expandedSvs <- merge(
            expandedSvs,
            unique(svRates[, .(group, Category1, Category2)]),
            all.x = TRUE,
            by = "group"
        )
    } else if(groupType == "bySample"){
        expandedSvs <- merge(
            expandedSvs,
            svRates[, .(uniqueId, Category1, Category2, Sample_ID)],
            all.x = TRUE,
            by.x = "group",
            by.y = "uniqueId"
        )
    } else {
        expandedSvs[, ":="(
            Category1 = 1,
            Category2 = 1
        )]
    }
    expandedSvs[order(Category1, Category2)]
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
    d <- svsByGroup[filter, .N, keyby = c("Category1", "Category2", "group", column)]
    d[, Frequency := N / sum(N), keyby = c("Category1", "Category2", "group")]
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
    groupType <- settings$get("Data", "Distribute")
    if(groupType != "allSamplesTogether"){
        plot$addLegend(
            legend = if(groupType == "byGroup") gsub("\n", ", ", groups) 
                     else getSampleNames(sampleUniqueIds = groups),
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
                abline(v = seq(-100, 100, 2),  col = "grey80")
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
            Project     = Project,
            Sample      = getSampleNames(sampleUniqueIds = uniqueId),   
            Group       = Category1Name,
            Category    = Category2Name,
            nSvs        = nSvs, 
            onTargetCoverage = onTargetCoverage,
            svRate      = svRate
        )]
    })
)

# ----------------------------------------------------------------------
# define bookmarking actions
# ----------------------------------------------------------------------
bookmarkObserver <- observe({
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
    bookmarkObserver$destroy()
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
