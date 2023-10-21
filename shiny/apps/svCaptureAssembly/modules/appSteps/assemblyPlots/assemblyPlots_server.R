#----------------------------------------------------------------------
# server components for the assemblyPlots appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
assemblyPlotsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module and settings
#----------------------------------------------------------------------
module <- 'assemblyPlots'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
denominatorColumn <- reactive({
    x <- "kbp_on_target"
    if(isTruthy(settings$get("Data","Enforce_Min_Family_Size"))) x <- paste(x, "filtered", sep = "_")
    if(isTruthy(settings$get("Data","Adjust_Read_Lengths")))     x <- paste(x, "effective", sep = "_")
    x
})
mergeDoses  <- reactive({ isTruthy(settings$get("Data","Merge_Doses")) })
mergeClones <- reactive({ isTruthy(settings$get("Data","Merge_Clones")) })
APC <- CONSTANTS$assemblyPlots
cacheCreateLevel <- "asNeeded"

#----------------------------------------------------------------------
# toggle data processing to allow manipulation of selections without too-frequence updates
#----------------------------------------------------------------------
observeEvent(input$suspendDataProcessing, {
    updateButton(
        session, 
        session$ns("suspendDataProcessing"), 
        label = if(input$suspendDataProcessing) "Allow Data Processing" else "Suspend Data Processing",
        style = if(input$suspendDataProcessing) "success" else "danger"
    )
})

#----------------------------------------------------------------------
# saving plot sets
#----------------------------------------------------------------------
workingId <- NULL # set to a plot id when editing a previously saved set
sendFeedback <- function(x, ...) output$savePlotSetFeedback <- renderText(x)
getPlotSetName <- function(id){
    name <- savedPlotSets$names[[id]] # user name overrides
    if(is.null(name)) savedPlotSets$list[[id]]$Name else name
}
getPlotSetNames <- function(rows = TRUE){
    sapply(names(savedPlotSets$list)[rows], getPlotSetName)
}
savedPlotSetsTemplate <- data.table(
    Remove      = character(),
    Name        = character(),
    Source      = character(),
    Group_By    = character(),
    Required    = character(),
    Prohibited  = character(),
    SV_Types    = character(),
    Projects    = character()
)
savedPlotSets <- summaryTableServer(
    id = 'savedPlotSets', # NOT ns(id) when nesting modules!
    parentId = id,
    stepNumber = options$stepNumber,
    stepLocks = locks[[id]],
    sendFeedback = sendFeedback,
    template = savedPlotSetsTemplate,
    type = 'shortList',
    remove = list(
        message = "Remove this set of saved plots?",
        name = getPlotSetName
    ),
    names = list(
        get = getPlotSetNames,
        source = id
    )
) 
observeEvent(input$savePlotSet, {
    sourceId <- sourceId()
    req(sourceId, input$savePlotSet, input$savePlotSet != 0)
    d <- list( # plot-defining metadata, shown on Saved Plots table; these define the data available to the plot
        Name = paste("Plot #", length(savedPlotSets$list) + 1),
        Source = getSourceFilePackageName(sourceId), # use the source name, not its unique ID, to allow sample additions to saved plots
        Group_By     = input$Group_By,
        Required     = input$Required,
        Prohibited   = input$Prohibited,
        SV_Types     = input$SV_Types,
        Projects     = input$Projects
    )
    r <- initializeRecordEdit(d, workingId, savedPlotSets$list, 'Plot Set', 'plot set', sendFeedback)
    d <- c(d, list( # non-definining formatting attributes saved with plot but not displayed on Saved Plots table
        conditions_svFrequencies = conditions_svFrequencies(), 
        groups_svFrequencies = groups_svFrequencies(),
        settings_svFrequencies = svFrequenciesPlot$settings$all_(),

        conditions_microhomology = conditions_microhomology(), 
        groups_microhomology = groups_microhomology(),
        settings_microhomology = microhomologyPlot$settings$all_(),

        conditions_endpoints = conditions_endpoints(), 
        groups_endpoints = groups_endpoints(),
        settings_endpoints = endpointsPlot$settings$all_(),

        conditions_svSizes = conditions_svSizes(), 
        groups_svSizes = groups_svSizes(),
        settings_svSizes = svSizesPlot$settings$all_()
    ))
    saveEditedRecord(d, workingId, savedPlotSets, r)
    workingId <<- NULL
})
addDataListObserver(module, savedPlotSetsTemplate, savedPlotSets, function(r, id){
    dt <- data.table(
        Remove = '', 
        Name   = '',
        Source = sub(".svCapture.assemble", "", r$Source)
    )
    for(x in c("Group_By","Required","Prohibited","SV_Types","Projects")) dt[[x]] <- paste(r[[x]], collapse = "<br>")
    dt
})
observeEvent(savedPlotSets$selected(), {
    plotSetI <- savedPlotSets$selected()
    clearLoadingPlotSet <- function(...){
        loadingPlotSet <<- NULL
        stopSpinner(session)
    } 
    abortSavedPlotSet <- function(...){
        sourceIdOverride(NA)
        workingId <<- NULL
        clearLoadingPlotSet()
    }  
    if(isTruthy(plotSetI)){
        plotSet <- savedPlotSets$list[[plotSetI]]
        sources <- app$upload$outcomes$sources()
        sourceI <- which(sapply(sources, function(source) source$manifest$Project == plotSet$Source))[1]
        if(isTruthy(sourceI)){
            startSpinner(session, message = "loading saved plots")
            sourceIdOverride(NA)
            setTimeout(function(....){ # let the prior plot clear first
                loadingPlotSet <<- plotSet
                waitFor(svFrequenciesPlotTrigger, clearLoadingPlotSet, delay = 100)

                svFrequenciesPlot$settings$replace(plotSet$settings_svFrequencies)
                microhomologyPlot$settings$replace(plotSet$settings_microhomology)
                endpointsPlot$settings$replace(plotSet$settings_endpoints)
                svSizesPlot$settings$replace(plotSet$settings_svSizes)

                sourceIdOverride(names(sources)[sourceI]) # thus, look up the current source of the saved name, which may be updated from the original save
                setTimeout(clearLoadingPlotSet, delay = 5000) # in case svFrequenciesPlotTrigger never fires due to failed plot
            }, delay = 100)
        } else abortSavedPlotSet()
    } else abortSavedPlotSet()
})

#----------------------------------------------------------------------
# assembly selection and loading
#----------------------------------------------------------------------
sourceIdOverride <- reactiveVal(NULL)
sourceId <- dataSourceTableServer(
    "dataSource", 
    selection = "single",
    sourceIdOverride = sourceIdOverride
)
assembly <- reactive({
    sourceId <- sourceId()
    denominatorColumn <- denominatorColumn()
    req(sourceId, denominatorColumn)
    rdsFile <- getSourceFilePath(sourceId, "assembly")
    req(file.exists(rdsFile))
    assemblyCache$get(
        'assembly', 
        permanent = TRUE,
        from = "ram",
        create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        keyObject = list(
            rdsFile = rdsFile,
            denominatorColumn = denominatorColumn,
            internalUseSampleColumns = internalUseSampleColumns
        ), 
        createFn = function(...) {
            startSpinner(session, message = "loading assembly")
            assembly <- readRDS(rdsFile)
            for(col in names(assembly$samples)){
                if(col %in% internalUseSampleColumns) next
                values <- assembly$samples[[col]] # ensure that all empty/zero-dose cells have the value "-"
                assembly$samples[[col]][is.na(values) | is.null(values) | values == "" | values == "0"] <- "-"

                ######################## TODO: fix this in samples list and rerun assembly, then delete
                assembly$samples[[col]][values == "bulk"] <- "-"

                values <- assembly$samples[[col]] # remove columns that never vary over all samples
                if(length(unique(values)) == 1) assembly$samples[[col]] <- NULL
            }
            assembly$samples$coverage <- as.integer(assembly$samples[[denominatorColumn]] / assembly$samples$sumTargetLens * 1000)
            stopSpinner(session)
            assembly
        }
    )$value
})

#----------------------------------------------------------------------
# selection of columns and values to group and plot, applies to all plots
#----------------------------------------------------------------------
groupableColumns <- reactive({
    cols <- names(assembly()$samples)
    cols[!(cols %in% internalUseSampleColumns)]
})
updateColumnSelectors <- function(choices, selected = NULL){
    default <- list(
        Group_By    = character(), 
        Required    = character(), 
        Prohibited  = character(),
        SV_Types    = "deletion"
    )
    if(is.null(selected)) selected <- default
    for(x in names(default)){
        choices_ <- switch(x,
            SV_Types = c("deletion", "duplication", "inversion", "translocation"),
            choices
        )
        updateCheckboxGroupInput(
            session = session,
            inputId = x,
            choices = choices_,
            selected = selected[[x]][selected[[x]] %in% choices_],
            inline = TRUE
        )
    } 
}
observeEvent(groupableColumns(), { updateColumnSelectors(groupableColumns(), loadingPlotSet) })

#----------------------------------------------------------------------
# optional table of all samples, collapsed on page load
#----------------------------------------------------------------------
allSamplesTable <- bufferedTableServer(
    "allSamples",
    id,
    input,
    reactive({ 
        samples <- assembly()$samples
        groupableColumns <- groupableColumns()
        req(samples, groupableColumns)
        samples[, .SD, .SDcols = c("project", "sample", groupableColumns, "coverage", svTypeColumns)] 
    }),
    selection = 'none',
    options = list(),
    filterable = TRUE
)

#----------------------------------------------------------------------
# working table of all samples matching grouping filters (not displayed to user)
#----------------------------------------------------------------------
groupedSamples <- reactive({
    samples <- assembly()$samples
    groupableColumns <- groupableColumns()
    req(samples, groupableColumns)
    groupedColumns <- input$Group_By
    if(length(groupedColumns) == 0) groupedColumns = groupableColumns
    ungroupedColumns <- groupableColumns[!(groupableColumns %in% groupedColumns)]
    assemblyCache$get(
        'groupedSamples', 
        permanent = TRUE,
        from = "ram",
        create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        keyObject = list(
            samples = samples,
            groupableColumns = groupableColumns,
            groupedColumns = groupedColumns,
            ungroupedColumns = ungroupedColumns,
            required = input$Required,
            prohibited = input$Prohibited,
            SV_Types = input$SV_Types
        ), 
        createFn = function(...) {
            startSpinner(session, message = "getting samples")
            for(column in unique(c(ungroupedColumns, input$Prohibited))){
                I <- samples[[column]] == "-" # all ungrouped/prohibited columns must be "-"
                samples <- samples[I]
            }
            for(column in input$Required){
                I <- samples[[column]] == "-" # all required columns must not be "-"
                samples <- samples[!I]
            }
            isVariable <- sapply(groupedColumns, function(column) length(unique(samples[[column]])) > 1) # only report informative columns
            variableColumns <- groupedColumns[isVariable]
            x <- samples[, .SD, .SDcols = c("project", "sample", variableColumns, "coverage", input$SV_Types)]
            if(mergeDoses())  for(column in variableColumns) x[[column]] <- dosesToLogical(x, column)
            if(mergeClones()) for(column in variableColumns) x[[column]] <- clonesToTargets(x, column)
            stopSpinner(session)
            x
        }
    )$value
})

# ----------------------------------------------------------------------
# selection of projects to group and plot; allows quick dropping of all samples in a project
# ----------------------------------------------------------------------
updateProjectSelector <- function(selected = NULL){
    groupedSamples <- groupedSamples()
    req(groupedSamples)
    projects <- unique(groupedSamples$project) # by default, all projects are selected
    updateCheckboxGroupInput(
        session  = session,
        inputId  = "Projects",
        choices  = projects,
        selected = if(is.null(selected)) projects else selected$Projects,
        inline = TRUE
    )
}
observeEvent(groupedSamples(), { updateProjectSelector(loadingPlotSet) })

# ----------------------------------------------------------------------
# table of all samples after applying grouping and project filters
# ----------------------------------------------------------------------
groupedProjectSamples <- reactive({ # not a slow step, not worth caching
    groupedSamples <- groupedSamples()
    req(groupedSamples)
    startSpinner(session, message = "getting project samples")
    x <- groupedSamples[project %in% input$Projects]
    stopSpinner(session)
    x
})
groupedProjectSamplesTable <- bufferedTableServer( # called Matching Samples on screen, starts collapsed
    "groupedProjectSamples",
    id,
    input,
    groupedProjectSamples,
    # editBoxes = list(),
    selection = 'none',
    options = list()
)

# ----------------------------------------------------------------------
# aggregate groups automatically determined from groupedProjectSamples
# ----------------------------------------------------------------------
groupingCols <- reactive({
    cols <- names(groupedProjectSamples())
    cols[!(cols %in% c("project", "sample", "coverage", denominatorColumns, svTypeColumns))]    
})
groups <- reactive({
    req(!input$suspendDataProcessing)
    gps <- groupedProjectSamples()
    groupingCols <- groupingCols()
    req(gps, input$SV_Types, groupingCols)
    assemblyCache$get(
        'groups', 
        permanent = TRUE,
        from = "ram",
        create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        keyObject = list(
            groupedProjectSamples = gps,
            groupingCols = groupingCols,
            SV_Types = input$SV_Types
        ), 
        createFn = function(...) {
            startSpinner(session, message = "aggregating sample groups")
            gps$nSvs <- 0
            gps$svFrequency <- 0
            for(svType in input$SV_Types) {
                gps$nSvs <- gps$nSvs + gps[[svType]]
                gps$svFrequency <- gps$svFrequency + gps[[svType]] / gps$coverage
            }
            groups <- gps[, .(
                nProjects = length(unique(project)),
                nSamples = .N,
                coverage = sum(coverage),
                nSvs = sum(nSvs),
                meanFrequency = round(mean(svFrequency), 4),
                sdFrequency = round(sd(svFrequency), 4),
                svFrequencies = list(svFrequency),
                projects = list(project)
            ), by = groupingCols]
            groups <- setGroupLabels(groups, groupingCols)
            stopSpinner(session)
            groups
        }
    )$value
})
groupsTable <- bufferedTableServer(
    "groups",
    id,
    input,
    reactive({ 
        groups <- groups() 
        cols <- names(groups)
        groups[, .SD, .SDcols = cols[!(cols %in% c("svFrequencies","projects","groupLabels"))]]
    }),
    # editBoxes = list(),
    selection = 'none',
    options = list(
        paging = FALSE,
        searching = FALSE
    )
)

#----------------------------------------------------------------------
# svFrequencies plot
#----------------------------------------------------------------------
svFrequenciesType <- "svFrequencies"
conditions_svFrequencies <- reactiveVal()
output$conditions_svFrequencies <- renderConditionsBucket(session, groupingCols, svFrequenciesType, conditions_svFrequencies)
observeEvent(input$conditions_svFrequenciesOnChange, { conditions_svFrequencies(input$conditions_svFrequenciesOnChange) })
#----------------------------------------------------------------------
groups_svFrequencies <- reactiveVal(NULL)
output$groups_svFrequencies <- renderGroupsBucket(session, groups, svFrequenciesType, groups_svFrequencies)
observeEvent(input$groups_svFrequenciesOnChange, { groups_svFrequencies(input$groups_svFrequenciesOnChange) })
#----------------------------------------------------------------------
svFrequenciesPlotTrigger <- reactiveVal(1)
svFrequenciesData <- reactive({ # not a slow step, not worth caching
    req(!input$suspendDataProcessing)
    conditions_ <- conditions_svFrequencies()
    groups_ <- groups_svFrequencies()
    req(conditions_, groups_)
    nConditions <- length(conditions_)
    nGroups     <- length(groups_)
    mar <- APC$svFrequencies$mar
    mar[1] <- mar[1] + nConditions
    list(
        Plot_Frame = getAssemblyPlotFrame(
            svFrequenciesPlot, 
            svFrequenciesPlot$settings$get("Groups","Group_Width_Inches") * nGroups, 
            APC$svFrequencies$insideHeight, 
            mar
        ),
        mar = mar,
        conditions_ = conditions_,
        groups_ = groups_,
        nConditions = nConditions,
        nGroups = nGroups
    )
})
svFrequenciesPlotFrame <- plotFrameReactive(svFrequenciesData) 
svFrequenciesPlot <- staticPlotBoxServer(
    "svFrequenciesPlot",
    settings = assemblyPlotSettings$svFrequencies, 
    size = "m",
    Plot_Frame = svFrequenciesPlotFrame,
    create = function() {
        req(!input$suspendDataProcessing)
        groups <- groups()
        groupingCols <- groupingCols()
        req(groups, groupingCols)
        d <- svFrequenciesData()
        req(d)
        startSpinner(session, message = "rendering frequency plot")
        conditionsI <- sapply(d$conditions_, function(x) which(groupingCols == x))
        groupsI     <- sapply(d$groups_,     function(x) which(groups$groupLabels == x))
        groups <- groups[groupsI]
        uniqueProjects <- unique(unlist(groups$projects))
        colors <- CONSTANTS$plotlyColors[1:length(uniqueProjects)]
        names(colors) = uniqueProjects
        maxY <- max(
            unlist(groups$svFrequencies), 
            groups[, meanFrequency + sdFrequency * 2],
            na.rm = TRUE
        ) * 1.05
        par(mar = d$mar)
        svFrequenciesPlot$initializeFrame(
            xlim = c(0.5, d$nGroups + 0.5),
            ylim = c(0, maxY),
            xlab = "",
            ylab = "SV Frequency", # SVs / Target Coverage
            yaxs = "i",
            xaxt = "n"
        )
        rect( # make the bar plot
            1:d$nGroups - 0.35, 
            0, 
            1:d$nGroups + 0.35, 
            groups[, meanFrequency], 
            lty = 1, 
            lwd = 1,
            col = "grey80"
        )
        for(i in 1:d$nGroups){ # overplot individual data points on the bar plot
            lines(rep(i, 2), groups[i, meanFrequency + sdFrequency * c(-2, 2)])
            svFrequencies <- unlist(groups[i, svFrequencies])
            projects <- unlist(groups[i, projects])
            svFrequenciesPlot$addPoints(
                x = jitter2(svFrequencies, i - 0.25, i + 0.25),
                y = svFrequencies,
                col = sapply(projects, function(x) colors[[x]])
            )
        }
        assemblyPlotConditionsGrid(groupingCols, groups, conditionsI)
        assemblyPlotTitle(svFrequenciesPlot, sourceId)
        isolate({ svFrequenciesPlotTrigger( svFrequenciesPlotTrigger() + 1 ) })
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# SVs found in the working samples, to populate distribution plots
#----------------------------------------------------------------------
matchingSvs <- reactive({
    req(!input$suspendDataProcessing)
    assembly <- assembly()
    gps <- groupedProjectSamples() 
    groupingCols <- groupingCols()   
    req(assembly, gps, groupingCols, input$SV_Types)
    assemblyCache$get(
        'matchingSvs', 
        permanent = TRUE,
        from = "ram",
        create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        keyObject = list(
            svs = assembly$svs,
            groupedProjectSamples = gps,
            groupingCols = groupingCols,
            SV_Types = input$SV_Types
        ), 
        createFn = function(...) {
            startSpinner(session, message = "collecting matching SVs")
            x <- gps[, {
                assembly$svs[
                    PROJECT == project & 
                    grepl(paste0(",", sample, ","), SAMPLES) & 
                    JXN_TYPE %in% svTypesToSymbols(input$SV_Types)]
            }, by = c("project", "sample", groupingCols)]
            x <- setGroupLabels(x, groupingCols)
            stopSpinner(session)
            x
        }
    )$value
})

#----------------------------------------------------------------------
# microhomology distribution plot
#----------------------------------------------------------------------
microhomologyType <- "microhomology"
conditions_microhomology <- reactiveVal(NULL)
output$conditions_microhomology <- renderConditionsBucket(session, groupingCols, microhomologyType, conditions_microhomology)
observeEvent(input$conditions_microhomologyOnChange, { conditions_microhomology(input$conditions_microhomologyOnChange) })
#----------------------------------------------------------------------
groups_microhomology <- reactiveVal(NULL)
output$groups_microhomology <- renderGroupsBucket(session, groups, microhomologyType, groups_microhomology)
observeEvent(input$groups_microhomologyOnChange, { groups_microhomology(input$groups_microhomologyOnChange) })
#----------------------------------------------------------------------
microhomologyData <- reactive({
    req(!input$suspendDataProcessing)
    svs <- matchingSvs()    
    groupingCols <- groupingCols()
    conditions_ <- conditions_microhomology()
    groups_ <- groups_microhomology()
    req(svs, groupingCols, conditions_, groups_)
    assemblyCache$get(
        'microhomologyData', 
        permanent = TRUE,
        from = "ram",
        create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        keyObject = list(
            svs = svs,
            groupingCols = groupingCols,
            conditions_ = conditions_,
            groups_ = groups_,
            constants = APC$microhomology,
            microhomologyPlot$settings$all()
        ), 
        createFn = function(...) {
            startSpinner(session, message = "counting microhomologies")
            svs <- svs[groupLabels %in% groups_ & N_SPLITS > 0] # omit gap-only junctions
            svs[, insertSize := -MICROHOM_LEN]
            titleSuffix <- NULL
            if(length(conditions_) < length(groupingCols)) {
                svs <- setGroupLabels(svs, groupingCols, conditions_) # regroup a second time if some conditions were omitted
                groups_ <- unique(svs$groupLabels)
                titleSuffix <- getDroppedGroupLabels(svs, groupingCols, conditions_)
            }
            dt <- dcastSvsByGroup(svs, as.formula("insertSize ~ groupLabels"), groups_)
            mar <- APC$microhomology$mar
            mar[3] <- mar[3] + length(groups_)
            x <- list(
                Plot_Frame = getAssemblyPlotFrame(
                    microhomologyPlot,
                    APC$microhomology$insideWidth, 
                    APC$microhomology$insideHeight, 
                    mar
                ),
                mar = mar,
                dt = dt,
                titleSuffix = titleSuffix,
                groupCounts = svs[, .N, by = .(groupLabels)]
            )
            stopSpinner(session)
            x
        }
    )$value
})
#----------------------------------------------------------------------
microhomologyPlotFrame <- plotFrameReactive(microhomologyData) 
microhomologyPlot <- staticPlotBoxServer(
    "microhomologyPlot",
    settings = assemblyPlotSettings$microhomology, 
    size = "m",
    Plot_Frame = microhomologyPlotFrame,
    create = function() {
        d <- microhomologyData()
        req(d, ncol(d$dt) > 1)
        startSpinner(session, message = "rendering microhomology plot")
        groupLabels_ <- names(d$dt)[2:ncol(d$dt)]
        maxY <- d$dt[, max(.SD, na.rm = TRUE), .SDcols = groupLabels_] * 1.05
        nGroups <- length(groupLabels_)
        par(mar = d$mar)
        xlim <- c(
            microhomologyPlot$settings$get("Limits","Min_X_Value"),
            microhomologyPlot$settings$get("Limits","Max_X_Value")
        )
        microhomologyPlot$initializeFrame(
            xlim = xlim, # TODO: expose as setting
            ylim = c(0, maxY),
            xlab = "Insert Size (bp)",
            ylab = "Frequency",
            yaxs = "i"
        )
        abline(v = c(seq(-50, 50, 5), -1, -2), col = "grey")
        abline(v = 0) 
        lwd <- 1.5
        assemblyPlotLines(microhomologyPlot, d$dt, lwd)
        assemblyPlotGroupsLegend(microhomologyPlot, xlim, maxY, groupLabels_, lwd, d$groupCounts)
        assemblyPlotTitle(microhomologyPlot, sourceId, d$titleSuffix)
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# endpoints distribution plot
#----------------------------------------------------------------------
endpointsType <- "endpoints"
conditions_endpoints <- reactiveVal(NULL)
output$conditions_endpoints <- renderConditionsBucket(session, groupingCols, endpointsType, conditions_endpoints)
observeEvent(input$conditions_endpointsOnChange, { conditions_endpoints(input$conditions_endpointsOnChange) })
#----------------------------------------------------------------------
groups_endpoints <- reactiveVal(NULL)
output$groups_endpoints <- renderGroupsBucket(session, groups, endpointsType, groups_endpoints)
observeEvent(input$groups_endpointsOnChange, { groups_endpoints(input$groups_endpointsOnChange) })
#----------------------------------------------------------------------
endpointsData <- reactive({
    req(!input$suspendDataProcessing)
    svs <- matchingSvs()
    groupingCols <- groupingCols()
    conditions_ <- conditions_endpoints()
    groups_ <- groups_endpoints()
    targets <- assembly()$targets
    binSize <- endpointsPlot$settings$get("Bins","Bin_Size")
    req(svs, groupingCols, conditions_, groups_, targets, binSize)
    assemblyCache$get(
        'endpointsData', 
        permanent = TRUE,
        from = "ram",
        # create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        create = 'once', # 'asNeeded', 'once', 'always'
        keyObject = list(
            svs = svs,
            groupingCols = groupingCols,
            conditions_ = conditions_,
            groups_ = groups_,
            targets = targets,
            binSize = binSize,
            constants = APC$endpoints,
            endpointsPlot$settings$all()
        ), 
        createFn = function(...) {
            startSpinner(session, message = "counting endpoints")
            targets[, chromI := {
                x <- sub("chr", "", chrom)
                x <- ifelse(x == "X", 99, x)
                x <- ifelse(x == "Y", 100, x)
                as.integer(x)
            }]
            targets <- targets[order(chromI, paddedStart)]
            nTargets <- nrow(targets)
            svs <- svs[groupLabels %in% groups_]
            titleSuffix <- NULL
            if(length(conditions_) < length(groupingCols)) {
                svs <- setGroupLabels(svs, groupingCols, conditions_) # regroup a second time if some conditions were omitted
                groups_ <- unique(svs$groupLabels)
                titleSuffix <- getDroppedGroupLabels(svs, groupingCols, conditions_)
            }
            svs <- svs[, .(endpointBin = as.integer(c(POS_1, POS_2) / binSize) * binSize), by = .(TARGET_REGION, groupLabels)]
            mar <- APC$endpoints$mar
            mar[3] <- mar[3] + length(groups_)
            regions <- lapply(1:nTargets, function(i){
                svs <- svs[TARGET_REGION == targets[i, name]]
                mar <- mar
                if(i != nTargets) mar[1] <- APC$untitledAxisMar
                if(i != 1) mar[3] <- APC$noMar
                dt <- dcastSvsByGroup(
                    svs[, .SD, .SDcols = c("endpointBin", "groupLabels")], 
                    as.formula("endpointBin ~ groupLabels"), 
                    groups_, 
                    step = binSize
                )
                list(
                    bins = dt, 
                    mar = mar,
                    regionCount = nrow(svs)
                )
            })
            Plot_Frame <- getAssemblyPlotFrame(
                endpointsPlot,
                APC$endpoints$insideWidth, 
                APC$endpoints$insideHeight * nTargets + APC$untitledAxisMar * (nTargets - 1) / APC$linesPerInch, 
                mar
            )
            x <- list(
                Plot_Frame = Plot_Frame,
                regions = regions,
                targets = targets,
                nTargets = nTargets,
                heights = sapply(1:nTargets, function(i) {
                    APC$endpoints$insideHeight + sum(regions[[i]]$mar[c(1, 3)]) / APC$linesPerInch
                }) / Plot_Frame$Height_Inches,
                titleSuffix = titleSuffix,
                groupCounts = svs[, .N, by = .(groupLabels)]
            )
            stopSpinner(session)
            x
        }
    )$value
})
#----------------------------------------------------------------------
endpointsPlotFrame <- plotFrameReactive(endpointsData) 
endpointsPlot <- staticPlotBoxServer(
    "endpointsPlot",
    settings = assemblyPlotSettings$endpoints, 
    size = "m",
    Plot_Frame = endpointsPlotFrame,
    create = function() {
        d <- endpointsData()
        req(d)
        startSpinner(session, message = "rendering endpoints plot")
        layout(matrix(1:d$nTargets, ncol = 1), heights = d$heights)
        for(i in 1:d$nTargets){
            region <- d$regions[[i]]
            geneStart <- d$targets[i, geneStart]
            geneEnd <- d$targets[i, geneEnd]

            groupLabels_ <- names(region$bins)[2:ncol(region$bins)]
            nGroups <- length(groupLabels_)            
            maxY <- region$bins[, max(.SD, na.rm = TRUE), .SDcols = groupLabels_] * 1.05
            par(mar = region$mar, cex = 1)
            xlim <- c(
                d$targets[i, paddedStart + 1],
                d$targets[i, paddedEnd]
            ) / 1e6
            ylim <- c(0, maxY)
            endpointsPlot$initializeFrame(
                xlim = xlim, 
                ylim = ylim,
                ylab = "Frequency",
                yaxs = "i",
                xlab = if(i == d$nTargets) "SV Endpoint Coordinate (Mbp)" else ""
            )
            rect(geneStart / 1e6, ylim[1] * 0.98, geneEnd / 1e6, ylim[2], col = "grey95", border = NA)
            abline(v = d$targets[i, c(regionStart, regionEnd) / 1e6], col = "grey40")
            # abline(v = d$targets[i, c(paddedStart, paddedEnd) / 1e6], col = "grey40")
            label <- paste(
                d$targets[i, paste(
                    chrom, 
                    paste(
                        name, 
                        if(is.na(geneStrand)) "" else paste0("(", geneStrand, ")")
                    ),
                    sep = ", "
                )],
                paste(region$regionCount, "SVs"), 
                sep = "\n"
            )
            text (xlim[1], maxY * 0.75, label, pos = 4, offset = 0)
            lwd <- 1
            assemblyPlotLines(endpointsPlot, region$bins, lwd, scale = 1e6)
            if(i == 1) {
                assemblyPlotGroupsLegend(endpointsPlot, xlim, maxY, groupLabels_, lwd, d$groupCounts)
                assemblyPlotTitle(endpointsPlot, sourceId, d$titleSuffix)
            }
        }
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# sizes distribution plot
#----------------------------------------------------------------------
svSizesType <- "svSizes"
conditions_svSizes <- reactiveVal(NULL)
output$conditions_svSizes <- renderConditionsBucket(session, groupingCols, svSizesType, conditions_svSizes)
observeEvent(input$conditions_svSizesOnChange, { conditions_svSizes(input$conditions_svSizesOnChange) })
#----------------------------------------------------------------------
groups_svSizes <- reactiveVal(NULL)
output$groups_svSizes <- renderGroupsBucket(session, groups, svSizesType, groups_svSizes)
observeEvent(input$groups_svSizesOnChange, { groups_svSizes(input$groups_svSizesOnChange) })
#----------------------------------------------------------------------
svSizesData <- reactive({
    req(!input$suspendDataProcessing)
    svs <- matchingSvs()
    groupingCols <- groupingCols()
    conditions_ <- conditions_svSizes()
    groups_ <- groups_svSizes()
    binsPerLog <- svSizesPlot$settings$get("Bins","Bin_Per_Log")
    req(svs, groupingCols, conditions_, groups_, binsPerLog)
    assemblyCache$get(
        'svSizesData', 
        permanent = TRUE,
        from = "ram",
        create = cacheCreateLevel, # 'asNeeded', 'once', 'always'
        keyObject = list(
            svs = svs,
            groupingCols = groupingCols,
            conditions_ = conditions_,
            groups_ = groups_,
            binsPerLog = binsPerLog,
            constants = APC$svSizes,
            svSizesPlot$settings$all()
        ), 
        createFn = function(...) {
            startSpinner(session, message = "counting sv sizes")
            svs <- svs[groupLabels %in% groups_]
            titleSuffix <- NULL
            if(length(conditions_) < length(groupingCols)) {
                svs <- setGroupLabels(svs, groupingCols, conditions_) # regroup a second time if some conditions were omitted
                groups_ <- unique(svs$groupLabels)
                titleSuffix <- getDroppedGroupLabels(svs, groupingCols, conditions_)
            }
            svs <- svs[, .(sizeBin = as.integer(log10(SV_SIZE) * binsPerLog) / binsPerLog), by = .(groupLabels)]
            dt <- dcastSvsByGroup(
                svs[, .SD, .SDcols = c("sizeBin", "groupLabels")], 
                as.formula("sizeBin ~ groupLabels"), 
                groups_, 
                step = 1 / binsPerLog
            )
            mar <- APC$svSizes$mar
            mar[3] <- mar[3] + length(groups_)            
            x <- list(
                Plot_Frame = getAssemblyPlotFrame(
                    svSizesPlot, 
                    APC$svSizes$insideWidth, 
                    APC$svSizes$insideHeight, 
                    mar
                ),
                mar = mar,
                dt = dt,
                titleSuffix = titleSuffix,
                groupCounts = svs[, .N, by = .(groupLabels)]
            )
            stopSpinner(session)
            x
        }
    )$value
})
#----------------------------------------------------------------------
svSizesPlotFrame <- plotFrameReactive(svSizesData) 
svSizesPlot <- staticPlotBoxServer(
    "svSizesPlot",
    settings = assemblyPlotSettings$svSizes, 
    size = "m",
    Plot_Frame = svSizesPlotFrame,
    create = function() {
        d <- svSizesData()
        req(d)
        startSpinner(session, message = "rendering svSizes plot")
        groupLabels_ <- names(d$dt)[2:ncol(d$dt)]
        maxY <- d$dt[, max(.SD, na.rm = TRUE), .SDcols = groupLabels_] * 1.05
        nGroups <- length(groupLabels_)
        par(mar = d$mar)
        xlim <- log10(c(
            svSizesPlot$settings$get("Bins","Min_SV_Size"),
            svSizesPlot$settings$get("Bins","Max_SV_Size")
        ))
        svSizesPlot$initializeFrame(
            xlim = xlim,
            ylim = c(0, maxY),
            xlab = "Log10 SV Size (bp)",
            ylab = "Frequency",
            yaxs = "i",
            title = ""
        )
        abline(v = 0:10, col = "grey")
        lwd <- 1
        assemblyPlotLines(svSizesPlot, d$dt, lwd)
        assemblyPlotGroupsLegend(svSizesPlot, xlim, maxY, groupLabels_, lwd, d$groupCounts)     
        assemblyPlotTitle(svSizesPlot, sourceId, d$titleSuffix)
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
selfDestruct <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    savedPlotSets$list  <- bm$outcomes$plotSets
    savedPlotSets$names <- bm$outcomes$plotSetNames
    selfDestruct$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(
        plotSets     = reactive(savedPlotSets$list),
        plotSetNames = reactive(savedPlotSets$names)
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
