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
linesPerInch <- 8.571429
fontSize <- 7
nullMar <- 0.5
stdLabelMar <- 4.1

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
        settings_endpoints = endpointsPlot$settings$all_()
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
loadingPlotSet <- NULL
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
    assembly
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
dosesToLogical <- function(x, col){ # collapse a set of different numerical doses (0.1, 0.2) to +/-
    I <- x[[col]] == "-"
    isNumeric <- all(!is.na(suppressWarnings(as.numeric(x[[col]][!I]))))
    if(isNumeric) x[[col]][!I] <- "+"
    x[[col]]
}
clonesToTargets <- function(x, col){ # collapse a set of clone replicates by stripping the clone name and retaining only the target
    sapply(x[[col]], function(value) strsplit(value, ":")[[1]][1]) # expects format TARGET:CLONE/REPLICATE
}
groupedSamples <- reactive({
    samples <- assembly()$samples
    groupableColumns <- groupableColumns()
    req(samples, groupableColumns)
    startSpinner(session, message = "getting samples")
    groupedColumns <- input$Group_By
    if(length(groupedColumns) == 0) groupedColumns = groupableColumns
    ungroupedColumns <- groupableColumns[!(groupableColumns %in% groupedColumns)]
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
    x
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
groupedProjectSamples <- reactive({
    groupedSamples <- groupedSamples()
    req(groupedSamples)
    startSpinner(session, message = "getting project samples")
    groupedSamples[project %in% input$Projects]
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
enDash <- function(x) {
    x[x == "-"] <- "\u2013"
    x
}
setGroupLabels <- function(x, groupingCols, labelCols = NULL){ # groupingCols not in labelCols are grouped and labeled as "all"
    workingCols <- if(is.null(labelCols)) groupingCols else groupingCols[groupingCols %in% labelCols]
    allCols <- if(is.null(labelCols)) character() else paste(groupingCols[!(groupingCols %in% labelCols)], "all", sep = " = ")
    x[, groupLabels := sapply(1:nrow(x), function(i) {
        paste(c(sapply(1:length(workingCols), function(j) {
            paste(workingCols[j], enDash(x[i][[workingCols[j]]]), sep = " = ")
        }), allCols), collapse = ", ")
    })]
    x
}
groups <- reactive({
    gps <- groupedProjectSamples()
    groupingCols <- groupingCols()
    req(gps, input$SV_Types, groupingCols)
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
    setGroupLabels(groups, groupingCols)
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
# svFrequencies plot element ordering
#----------------------------------------------------------------------
getBucketList <- function(rankListId, labels){
    keepRankListId <- paste0(rankListId, "Keep")
    dropRankListId <- paste0(rankListId, "Drop")
    rankListOnChangeId <- session$ns(paste0(rankListId, "OnChange"))
    if(is.null(loadingPlotSet)){
        keepLabels <- labels
        dropLabels <-character()
    } else {
        x <- loadingPlotSet[[rankListId]]
        keepLabels <- if(is.null(x)) labels else x[x %in% labels]
        dropLabels <- labels[!(labels %in% keepLabels)]
    }
    list(
        order = keepLabels,
        ui = bucket_list(
            "",
            add_rank_list(
                text = "Show/Use/Group on Plot",
                labels = keepLabels,
                input_id = keepRankListId,
                css_id = keepRankListId,
                options = sortable_options(
                    group = rankListId,
                    multiDrag = TRUE,
                    onSort = sortable_js_capture_input(input_id = rankListOnChangeId)
                ),
                class = "default-sortable" 
            ),
            add_rank_list(
                text = "Omit from / merge on Plot",
                labels = dropLabels,
                input_id = dropRankListId,
                css_id = dropRankListId,
                options = sortable_options(
                    group = rankListId,
                    multiDrag = TRUE
                ),
                class = "default-sortable" 
            )
        ) 
    )
}
#----------------------------------------------------------------------
conditions_svFrequencies <- reactiveVal()
output$conditions_svFrequencies <- renderUI({
    groupingCols <- groupingCols()
    req(groupingCols)
    x <- getBucketList("conditions_svFrequencies", groupingCols)
    conditions_svFrequencies(x$order)
    x$ui
})
observeEvent(input$conditions_svFrequenciesOnChange, { conditions_svFrequencies(input$conditions_svFrequenciesOnChange) })
#----------------------------------------------------------------------
groups_svFrequencies <- reactiveVal(NULL)
output$groups_svFrequencies <- renderUI({
    groups <- groups()
    req(groups)
    x <- getBucketList("groups_svFrequencies", groups$groupLabels)
    groups_svFrequencies(x$order)
    x$ui
})
observeEvent(input$groups_svFrequenciesOnChange, { groups_svFrequencies(input$groups_svFrequenciesOnChange) })

#----------------------------------------------------------------------
# svFrequencies plot
#----------------------------------------------------------------------
jitter2 <- function(v, min, max){
    N <- length(v)
    width <- max - min
    if(N == 1) return(min + width / 2)
    if(N <= 4) {
        min <- min + width / 4
        max <- max - width / 4
    } 
    x <- seq(min, max, length.out = N)
    x[sample(N)]        
}
svFrequenciesPlotTrigger <- reactiveVal(1)
svFrequenciesDim <- list(
    insideHeight = 1,
    groupWidth = 0.3,
    mar = c(1.1, 8.1, nullMar, nullMar)
)
svFrequenciesData <- reactive({
    conditions_ <- conditions_svFrequencies()
    groups_ <- groups_svFrequencies()
    req(conditions_, groups_)
    nConditions <- length(conditions_)
    nGroups     <- length(groups_)
    mar <- svFrequenciesDim$mar
    mar[1] <- mar[1] + nConditions
    maiHorizonatal <- sum(mar[c(2, 4)]) / linesPerInch
    maiVertical    <- sum(mar[c(1, 3)]) / linesPerInch
    list(
        Plot_Frame = list(
            Width_Inches  = svFrequenciesPlot$settings$get("Width","Group_Width_Inches") * nGroups + maiHorizonatal, 
            Height_Inches = svFrequenciesDim$insideHeight + maiVertical,
            Font_Size = fontSize
        ),
        mar = mar,
        conditions_ = conditions_,
        groups_ = groups_,
        nConditions = nConditions,
        nGroups = nGroups
    )
})
svFrequenciesPlotFrame <- reactive({ 
    tryCatch({
        svFrequenciesData()$Plot_Frame
    }, error = function(e) list(
        Width_Inches  = 3, 
        Height_Inches = 3,
        Font_Size = fontSize
    ))
})
svFrequenciesPlot <- staticPlotBoxServer(
    "svFrequenciesPlot",
    # title = TRUE,
    settings = list(
        Width = list(
            Group_Width_Inches = list(
                type = "numericInput",
                value = svFrequenciesDim$groupWidth,
                min = 0.1,
                max = 1,
                step = 0.05
            )
        )
    ), 
    size = "m",
    Plot_Frame = svFrequenciesPlotFrame,
    create = function() {
        groups <- groups()
        groupingCols <- groupingCols()
        req(groups, groupingCols)
        d <- svFrequenciesData()
        req(d)
        startSpinner(session, message = "rendering frequency plot")
        d$conditions_ <- sapply(d$conditions_, function(x) which(groupingCols == x))
        d$groups_     <- sapply(d$groups_,     function(x) which(groups$groupLabels == x))
        groups <- groups[d$groups_]
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
        mtext( # labels for the condition grid rows
            gsub("_", " ", groupingCols[d$conditions_]), 
            side = 1, 
            line = 1:d$nConditions, 
            at = 0,
            adj = 1,
            cex = 0.9
        )
        for(i in 1:d$nConditions){ # fill the condition grid with values by row
            j <- d$conditions_[i]
            mtext(
                enDash(unlist(groups[, .SD, .SDcols = groupingCols[j]])), 
                side = 1, 
                line = i, 
                at = 1:d$nGroups
            )
        }
        isolate({ svFrequenciesPlotTrigger( svFrequenciesPlotTrigger() + 1 ) })
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# SVs found in the working samples, to populate distribution plots
#----------------------------------------------------------------------
svTypesToSymbols <- function(types){
    sapply(types, switch, deletion = "L", duplication = "D", inversion = "I", translocation = "T")
}
matchingSvs <- reactive({
    assembly <- assembly()
    gps <- groupedProjectSamples()    
    req(assembly, gps)
    startSpinner(session, message = "collecting matching SVs")
    x <- gps[, {
        assembly$svs[
            PROJECT == project & 
            grepl(paste0(",", sample, ","), SAMPLES) & 
            JXN_TYPE %in% svTypesToSymbols(input$SV_Types)]
    }, by = c("project", "sample", groupingCols())]
    setGroupLabels(x, groupingCols())
})

#----------------------------------------------------------------------
# microhomology distribution plot
#----------------------------------------------------------------------
conditions_microhomology <- reactiveVal(NULL)
output$conditions_microhomology <- renderUI({
    groupingCols <- groupingCols()
    req(groupingCols)
    x <- getBucketList("conditions_microhomology", groupingCols)
    conditions_microhomology(x$order)
    x$ui
})
observeEvent(input$conditions_microhomologyOnChange, { conditions_microhomology(input$conditions_microhomologyOnChange) })
#----------------------------------------------------------------------
groups_microhomology <- reactiveVal(NULL)
output$groups_microhomology <- renderUI({
    groups <- groups()
    req(groups)
    x <- getBucketList("groups_microhomology", groups$groupLabels)
    groups_microhomology(x$order)
    x$ui
})
observeEvent(input$groups_microhomologyOnChange, { groups_microhomology(input$groups_microhomologyOnChange) })
#----------------------------------------------------------------------
microhomologyDim <- list(
    insideWidth  = 1.5,
    insideHeight = 1,
    mar = c(stdLabelMar, stdLabelMar, 1.1, nullMar)
)
microhomologyData <- reactive({
    groupingCols <- groupingCols()
    svs <- matchingSvs()
    conditions_ <- conditions_microhomology()
    groups_ <- groups_microhomology()
    req(svs, conditions_, groups_)
    startSpinner(session, message = "counting microhomologies")
    svs <- svs[groupLabels %in% groups_ & N_SPLITS > 0] # omit gap-only junctions
    if(length(conditions_) < length(groupingCols)) {
        svs <- setGroupLabels(svs, groupingCols, conditions_) # regroup a second time if some conditions were omitted
        groups_ <- unique(svs$groupLabels)
    }
    dt <- dcast(
        svs,
        MICROHOM_LEN ~ groupLabels, 
        value.var = "MICROHOM_LEN", 
        fun.aggregate = length, 
        fill = 0
    )
    for(gl in groups_) dt[[gl]] <- dt[[gl]] / sum(dt[[gl]])
    mar <- microhomologyDim$mar
    mar[3] <- mar[3] + length(groups_)
    maiHorizonatal <- sum(mar[c(2, 4)]) / linesPerInch
    maiVertical    <- sum(mar[c(1, 3)]) / linesPerInch
    list(
        Plot_Frame = list(
            Width_Inches  = microhomologyDim$insideWidth + maiHorizonatal, 
            Height_Inches = microhomologyDim$insideHeight + maiVertical,
            Font_Size = fontSize
        ),
        mar = mar,
        dt = dt
    )
})
#----------------------------------------------------------------------
microhomologyPlotFrame <- reactive({ 
    tryCatch({
        microhomologyData()$Plot_Frame
    }, error = function(e) list(
        Width_Inches  = 3, 
        Height_Inches = 3,
        Font_Size = fontSize
    ))
})
microhomologyPlot <- staticPlotBoxServer(
    "microhomologyPlot",
    # title   = TRUE,
    settings = list(
        Limits = list(
            Min_X_Value = list(
                type = "numericInput",
                value = -10,
                min = -50, 
                max = 0,
                step = 5
            ),
            Max_X_Value = list(
                type = "numericInput",
                value = 15,
                min = 0, 
                max = 50,
                step = 5
            )
        )
    ), 
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
        for(i in 1:nGroups){ # overplot individual data points on the bar plot
            microhomologyPlot$addLines(
                x = -d$dt$MICROHOM_LEN, # thus, plotting inserts as positive values, microhomology as negative
                y = d$dt[[groupLabels_[i]]],
                col = CONSTANTS$plotlyColors[[i]],
                lwd = lwd
            )
        }
        microhomologyPlot$addMarginLegend(
            x = mean(xlim),
            xjust = 0.5,
            y = maxY,
            yjust = 0,
            legend = groupLabels_,
            col = unlist(CONSTANTS$plotlyColors[1:nGroups]),
            lty = 1,
            lwd = lwd,
            bty = "n",
            cex = 0.8
        )
        stopSpinner(session)
    }
)

#----------------------------------------------------------------------
# endpoints distribution plot
#----------------------------------------------------------------------
conditions_endpoints <- reactiveVal(NULL)
output$conditions_endpoints <- renderUI({
    groupingCols <- groupingCols()
    req(groupingCols)
    x <- getBucketList("conditions_endpoints", groupingCols)
    conditions_endpoints(x$order)
    x$ui
})
observeEvent(input$conditions_endpointsOnChange, { conditions_endpoints(input$conditions_endpointsOnChange) })
#----------------------------------------------------------------------
groups_endpoints <- reactiveVal(NULL)
output$groups_endpoints <- renderUI({
    groups <- groups()
    req(groups)
    x <- getBucketList("groups_endpoints", groups$groupLabels)
    groups_endpoints(x$order)
    x$ui
})
observeEvent(input$groups_endpointsOnChange, { groups_endpoints(input$groups_endpointsOnChange) })
#----------------------------------------------------------------------
endpointsDim <- list(
    insideWidth  = 3,
    insideHeight = 1,
    mar = c(stdLabelMar, stdLabelMar, 1.1, nullMar)
)
endpointsData <- reactive({
    groupingCols <- groupingCols()
    svs <- matchingSvs()
    conditions_ <- conditions_endpoints()
    groups_ <- groups_endpoints()
    targets <- assembly()$targets
    nTargets <- nrow(targets)
    binSize <- endpointsPlot$settings$get("Bins","Bin_Size")
    req(svs, conditions_, groups_, targets, binSize)
    startSpinner(session, message = "counting endpoints")
    svs <- svs[groupLabels %in% groups_]
    if(length(conditions_) < length(groupingCols)) {
        svs <- setGroupLabels(svs, groupingCols, conditions_) # regroup a second time if some conditions were omitted
        groups_ <- unique(svs$groupLabels)
    }
    svs <- svs[, .(endpointBin = as.integer(c(POS_1, POS_2) / binSize) * binSize), by = .(TARGET_REGION, groupLabels)]
    mar <- endpointsDim$mar
    mar[3] <- mar[3] + length(groups_)
    maiHorizonatal <- sum(mar[c(2, 4)]) / linesPerInch
    maiVertical    <- sum(mar[c(1, 3)]) / linesPerInch
    regions <- lapply(1:nTargets, function(i){
        mar <- mar
        if(i != 1) mar[3] <- 0.1
        if(i != nTargets) mar[1] <- 0.1
        list(
            bins = dcast(
                svs[TARGET_REGION == targets[i, name]],
                endpointBin ~ groupLabels, 
                value.var = "endpointBin", 
                fun.aggregate = length, 
                fill = 0
            ),
            mar = mar
        )
        # for(gl in groups_) dt[[gl]] <- dt[[gl]] / sum(dt[[gl]])            
    })
    list(
        Plot_Frame = list(
            Width_Inches  = endpointsDim$insideWidth  * nTargets + maiHorizonatal, 
            Height_Inches = endpointsDim$insideHeight * nTargets + maiVertical,
            Font_Size = fontSize
        ),
        regions = regions,
        targets = targets,
        nTargets = nTargets
    )
#  $ TARGET_REGION             : chr  "WWOX" "WWOX" "WWOX" "WWOX" ...
#  $ SV_SIZE                   : num  1501 6331 1558 3609 1217 ...
#  $ POS_1                     : int  78494093 78494096 78498346 78498233 78498244 78498393 78499200 78498321 78497952 78496939 ...
#  $ POS_2 

# Classes ‘data.table’ and 'data.frame':  3 obs. of  10 variables:
#  $ genome        : chr  "hg38" "hg38" "hg38"
#  $ targetsBed    : chr  "/nfs/turbo/umms-glover/data/svCapture/resources/capture
# _targets/HF1_capture_targets.hg38.bed" "/nfs/turbo/umms-glover/data/svCapture/re
# sources/capture_targets/HF1_capture_targets.hg38.bed" "/nfs/turbo/umms-glover/da
# ta/svCapture/resources/capture_targets/HF1_capture_targets.hg38.bed"
#  $ chrom         : chr  "chr1" "chr10" "chr7"
#  $ regionStart   : int  71684317 51540240 78460683
#  $ regionEnd     : int  71934317 51790240 78710684
#  $ name          : chr  "NEGR1" "PRKG1" "MAGI2"
#  $ regionPadding : int  800000 800000 800000
#  $ paddedStart   : int  70884317 50740240 77660683
#  $ paddedEnd     : int  72734317 52590240 79510684
#  $ paddedSequence: chr  "GAAGCTGGGAATATCATGGGTAATAACTCCCATGATTATGTTACATTATATTGCC
# AAAGGTATTTTCCAGATGCAATTAAGTTTCCTAGTCAACTGATTTTGAGTTAATTG"| __truncated__ "AGGCGA
# CTCGGACCTGAGCCGAATCTGTTGACCCCAAATTGTGCTTTTCCCACCAAGAGCAAAAGAAAGAGAAACATTAGTACAGC
# TTCGGAACTAAAATATAGTAGAGAA"| __truncated__ "GCAGTGGCTCACATCTATAATCCCAGCACTTTGGGAG
# GGCAAAGCAGGCAGATCGTGTGAGCTCAGGAGTTCGAAACCAGCCTGGCCAACGTGGTGAAACTCTGTCTCTGC"| __t
# runcated__
})
#----------------------------------------------------------------------
endpointsPlotFrame <- reactive({ 
    tryCatch({
        endpointsData()$Plot_Frame
    }, error = function(e) list(
        Width_Inches  = 3, 
        Height_Inches = 3,
        Font_Size = fontSize
    ))
})
endpointsPlot <- staticPlotBoxServer(
    "endpointsPlot",
    # title   = TRUE,
    settings = list(
        Bins = list(
            Bin_Size = list(
                type = "numericInput",
                value = 1000,
                min = 100, 
                max = 100000,
                step = 100
            )
        )
    ), 
    size = "m",
    Plot_Frame = endpointsPlotFrame,
    create = function() {
        d <- endpointsData()
        req(d)
        startSpinner(session, message = "rendering endpoints plot")
        layout(matrix(1:d$nTargets, ncol = 1)) # heights 
        for(i in 1:d$nTargets){
            region <- d$regions[[i]]

            dstr(region)

            groupLabels_ <- names(region$bins)[2:ncol(region$bins)]
            maxY <- region$bins[, max(.SD, na.rm = TRUE), .SDcols = groupLabels_] * 1.05
            nGroups <- length(groupLabels_)
            par(mar = region$mar, cex = 1)
            xlim <- c(
                d$targets[i, paddedStart + 1],
                d$targets[i, paddedEnd]
            )
            endpointsPlot$initializeFrame(
                xlim = xlim, # TODO: expose as setting
                ylim = c(0, maxY),
                xlab = "Coordinate (bp)",
                ylab = "Frequency",
                yaxs = "i"
            )
            # abline(v = c(seq(-50, 50, 5), -1, -2), col = "grey")
            # abline(v = 0) 
            lwd <- 1.5
            # for(i in 1:nGroups){ # overplot individual data points on the bar plot
            #     endpointsPlot$addLines(
            #         x = -d$dt$MICROHOM_LEN, # thus, plotting inserts as positive values, microhomology as negative
            #         y = d$dt[[groupLabels_[i]]],
            #         col = CONSTANTS$plotlyColors[[i]],
            #         lwd = lwd
            #     )
            # }
            # endpointsPlot$addMarginLegend(
            #     x = mean(xlim),
            #     xjust = 0.5,
            #     y = maxY,
            #     yjust = 0,
            #     legend = groupLabels_,
            #     col = unlist(CONSTANTS$plotlyColors[1:nGroups]),
            #     lty = 1,
            #     lwd = lwd,
            #     bty = "n",
            #     cex = 0.8
            # )            
        }
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
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
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
