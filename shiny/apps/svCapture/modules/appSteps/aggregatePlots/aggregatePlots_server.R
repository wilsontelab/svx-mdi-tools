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
        file.path(app$sources$suiteGlobalDir, "sv_filters.yml"),
        id
    ),
    fade = FALSE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues() # logical failure vectors keyed as [[sampleSet]]

#----------------------------------------------------------------------
# generate the list of all filtered SVs from all selected samples
#----------------------------------------------------------------------
targetClasses <- reactive({ SVX$targetClasses[[settings$SV_Filters()$Target_Class$value]] })
workingSvs <- reactive({
    getWorkingSvs(settings, sampleSelector, parseSamples = TRUE)
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
    ps <- workingSvs()[, PROJECT_SAMPLES]
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
        get <- svRatesPlot$settings$get
        svRates <- svRates()[order(rev(1:.N))]
        req(svRates)

        svRates[, groupFactor := factor(group, unique(group))]
        # dstr(svRates)

        # barplot(as.matrix(x), beside = TRUE, horiz = TRUE)
        par(mar = c(
            get("Plot_Frame", "Bottom_Margin"), 
            get("Plot_Frame", "Left_Margin"), 
            get("Plot_Frame", "Top_Margin"), 
            get("Plot_Frame", "Right_Margin")
        ))
        xi <- as.integer(svRates$groupFactor)
        # plot(
        #     jitter(xi, amount = 0.1),
        #     svRates$svRate,
        #     typ = "p",

        #     xlim = c(0, max(as.integer(svRates$groupFactor))) + 0.5,
        #     xlab = "",
        #     ylim = c(0, max(svRates$svRate) * 1.1),
        #     ylab = "SV Rate",
        #     xaxt = "n"
        # )
        plot(
            NA, NA, 
            typ = "p",
            main = get("Plot_Frame", "Title"),

            ylim = c(0, max(as.integer(svRates$groupFactor))) + 0.5,
            ylab = "",
            xlim = c(0, max(svRates$svRate) * 1.1),
            xlab = "SV Rate",
            yaxt = "n"
        )
        abline(h = c(0, xi) + 0.5, col = "grey80")

        # dstr(CONSTANTS$plotlyColors)
        # dstr(as.integer(factor(svRates$replicate)))
        # dstr(CONSTANTS$plotlyColors[as.integer(factor(svRates$replicate))])
        points(
            svRates$svRate,            
            jitter(xi, amount = 0.1),
            pch = get("Points_and_Lines", "Point_Type"),
            cex = get("Points_and_Lines", "Point_Size"),
            col = unlist(CONSTANTS$plotlyColors[as.integer(factor(svRates$replicate))])
        )
        axis(
            2, 
            at = xi, 
            labels = svRates$group,
            las = 1
            # , 
            # tick = TRUE, 
            # line = NA,
            # pos = NA, 
            # outer = FALSE, 
            # font = NA, 
            # lty = "solid",
            # lwd = 1, 
            # lwd.ticks = lwd, 
            # col = NULL, 
            # col.ticks = NULL,
            # hadj = NA, 
            # padj = NA, 
            # gap.axis = NA
        )

        # barplot(svRate ~ group, assignments, beside = TRUE, horiz = TRUE)

        #     xaxs = "i", 
        #     yaxs = "i",
        #     xaxt = "n"
    }
)

#----------------------------------------------------------------------
# Microhomology Length Distribution plot
#----------------------------------------------------------------------
svsByGroup <- reactive({
    maxSamples <- settings$SV_Filters()$Max_Samples_With_SV$value
    req(maxSamples == 1)
    svRates <- svRates()
    req(svRates)
    setkey(svRates, uniqueId)        
    svs <- workingSvs()[, .(
        plotted = JXN_BASES != "*", # MICROHOM_LEN meaningless if not a sequenced junction
        MICROHOM_LEN, 
        SV_SIZE,
        group = svRates[unlist(PROJECT_SAMPLES), group]
    )]    
})
plotFrequencyDistribution <- function(
    column, xlab, xlim, 
    modifyData = function(...) NULL,  
    modifyPlotBase = function(...) NULL
){
    svsByGroup <- svsByGroup()
    req(svsByGroup)
    modifyData(svsByGroup)
    d <- svsByGroup[plotted == TRUE, .N, keyby = c("group", column)]
    d[, Frequency := N / sum(N), keyby = group]
    plot(
        NA, NA, typ = "n",
        xlim = xlim,
        ylim = c(0, max(d$Frequency) * 1.1),
        xlab = xlab,
        ylab = "Frequency"
    )
    modifyPlotBase(d)
    groups <- unique(d$group)
    for(i in seq_along(groups)){
        dd <- d[group == groups[i]]
        lines(dd[[column]], dd$Frequency, col = CONSTANTS$plotlyColors[[i]])
    }        
}
sizeDistribution <- staticPlotBoxServer(
    'sizeDistribution',
    legend = TRUE,
    immediate = TRUE,
    create = function(...){
        plotFrequencyDistribution(
            column = "SV_SIZE", 
            xlab = "log10 SV Size (bp)",
            xlim = log10(c(1000, 1e6)),
            modifyData = function(svsByGroup) svsByGroup[, SV_SIZE := round(log10(SV_SIZE), 1)]
            # ,
            # modifyPlotBase = function(...){
            #     abline(v = seq(-100, 100, 10), col = "grey60")
            #     abline(v = 0)
            # }
        )
    }
)
microhomologyDistribution <- staticPlotBoxServer(
    'microhomologyDistribution',
    legend = TRUE,
    immediate = TRUE,
    create = function(...){
        plotFrequencyDistribution(
            column = "MICROHOM_LEN", 
            xlab = "Microhomology Length (bp)",
            xlim = c(-20, 20),
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
            Group = Category1Name,
            Category = Category2Name,
            Project = Project,
            Sample = Sample_ID,
            nSvs = nSvs, 
            onTargetCoverage = onTargetCoverage,
            svRate = svRate
        )]
    })
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
        svRatesPlot$settings$replace(bm$outcomes$svRatesPlotPlotSettings)
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
        svRatesPlotPlotSettings  = svRatesPlot$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
