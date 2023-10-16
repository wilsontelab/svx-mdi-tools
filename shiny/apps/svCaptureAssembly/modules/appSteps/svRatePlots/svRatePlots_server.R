#----------------------------------------------------------------------
# server components for the svRatePlots appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svRatePlotsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module and settings
#----------------------------------------------------------------------
module <- 'svRatePlots'
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
mergeDoses <- reactive({
    isTruthy(settings$get("Data","Merge_Doses"))
})
mergeClones <- reactive({
    isTruthy(settings$get("Data","Merge_Clones"))
})

#----------------------------------------------------------------------
# saving figures
#----------------------------------------------------------------------
workingId <- NULL # set to a plot id when editing a previously saved set
sendFeedback <- function(x, ...){
    output$savePlotFeedback <- renderText(x)
}
getPlotName <- function(id){
    name <- savedPlots$names[[id]] # user name overrides
    if(is.null(name)) savedPlots$list[[id]]$Name else name
}
getPlotNames <- function(rows = TRUE){
    sapply(names(savedPlots$list)[rows], getPlotName)
}
savedPlotsTemplate <- data.table(
    Remove      = character(),
    Name        = character(),
    Source      = character(),
    Group_By    = character(),
    Required    = character(),
    Prohibited  = character(),
    SV_Types    = character(),
    Projects    = character()
)
savedPlots <- summaryTableServer(
    id = 'savedPlots', # NOT ns(id) when nesting modules!
    parentId = id,
    stepNumber = options$stepNumber,
    stepLocks = locks[[id]],
    sendFeedback = sendFeedback,
    template = savedPlotsTemplate,
    type = 'shortList',
    remove = list(
        message = "Remove this saved plot?",
        name = getPlotName
    ),
    names = list(
        get = getPlotNames,
        source = id
    )
) 
observeEvent(input$savePlot, {
    sourceId <- sourceId()
    req(
        sourceId, 
        input$savePlot, 
        input$savePlot != 0
    )   
    d <- list( # plot-defining metadata, shown on Saved Plots table; these define the data available to the plot
        Name = paste("Plot #", length(savedPlots$list) + 1),
        Source = getSourceFilePackageName(sourceId), # use the source name, not its unique ID, to allow sample additions to saved plots
        Group_By     = input$Group_By,
        Required     = input$Required,
        Prohibited   = input$Prohibited,
        SV_Types     = input$SV_Types,
        Projects     = input$Projects
    )
    r <- initializeRecordEdit(d, workingId, savedPlots$list, 'Plot', 'plot', sendFeedback)
    d <- c(d, list( # non-definining formatting attributes saved with plot but not displayed on Saved Plots table
        width = 3, # TODO: get from current settings
        conditionOrder = conditionOrder(), 
        groupOrder = groupOrder()
    ))
    saveEditedRecord(d, workingId, savedPlots, r)
    workingId <<- NULL
})
addDataListObserver(module, savedPlotsTemplate, savedPlots, function(r, id){
    dt <- data.table(
        Remove = '', 
        Name   = '',
        Source = sub(".svCapture.assemble", "", r$Source)
    )
    for(x in c("Group_By","Required","Prohibited","SV_Types","Projects")) dt[[x]] <- paste(r[[x]], collapse = "<br>")
    dt
})
loadingPlot <- NULL
observeEvent(savedPlots$selected(), {
    plotI <- savedPlots$selected()
    if(isTruthy(plotI)){
        plot <- savedPlots$list[[plotI]]
        sources <- app$upload$outcomes$sources()
        sourceI <- which(sapply(sources, function(source) source$manifest$Project == plot$Source))[1]
        req(sourceI)

        loadingPlot <<- plot
        startSpinner(session, message = "loading saved plot")
        sourceIdOverride(names(sources)[sourceI]) # thus, look up the current source of the saved name, which may be updated from the original save

        # change to waitFor of initial plotting in addition to timeout
        setTimeout(function(...){
            loadingPlot <<- NULL
            stopSpinner(session)
        }, delay = 1000)

    } else {
        sourceIdOverride(NA)
        workingId <<- NULL
    }
})

#----------------------------------------------------------------------
# assembly selection
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

        ########################
        assembly$samples[[col]][values == "bulk"] <- "-"

        values <- assembly$samples[[col]] # remove columns that never vary over all samples
        if(length(unique(values)) == 1) assembly$samples[[col]] <- NULL
    }
    assembly$samples$coverage <- as.integer(assembly$samples[[denominatorColumn]] / assembly$samples$sumTargetLens * 1000)
    assembly
})

#----------------------------------------------------------------------
# selection of columns and values to group and plot
#----------------------------------------------------------------------
# groupableColumnsTrigger <- reactiveVal(1)
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
observeEvent(groupableColumns(), {
    updateColumnSelectors(groupableColumns(), loadingPlot)
    # groupableColumnsTrigger( groupableColumnsTrigger() + 1 )
})

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
dosesToLogical <- function(x, col){ # convert a set of different numerical doses (0.1,0.2) 
    I <- x[[col]] == "-"
    isNumeric <- all(!is.na(suppressWarnings(as.numeric(x[[col]][!I]))))
    if(isNumeric) x[[col]][!I] <- "+"
    x[[col]]
}
clonesToTargets <- function(x, col){
    sapply(x[[col]], function(value) strsplit(value, ":")[[1]][1])
}
# groupedSamplesTrigger <- reactiveVal(1)
groupedSamples <- reactive({
    samples <- assembly()$samples
    groupableColumns <- groupableColumns()
    req(samples, groupableColumns)
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
    # isolate({ groupedSamplesTrigger( groupedSamplesTrigger() + 1 ) })
    x
})

# ----------------------------------------------------------------------
# selection of projects to group and plot
# ----------------------------------------------------------------------
# projectsTrigger <- reactiveVal(1)
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
observeEvent(groupedSamples(), {
    updateProjectSelector(loadingPlot)
    # projectsTrigger( projectsTrigger() + 1 )
})

# ----------------------------------------------------------------------
# table of all samples after applying grouping and project filters
# ----------------------------------------------------------------------
groupedProjectSamples <- reactive({
    groupedSamples <- groupedSamples()
    req(groupedSamples)
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
# groups automatically determined from groupedProjectSamples
# ----------------------------------------------------------------------
groupingCols <- reactive({
    cols <- names(groupedProjectSamples())
    cols[!(cols %in% c("project", "sample", "coverage", denominatorColumns, svTypeColumns))]    
})
groups <- reactive({
    gps <- groupedProjectSamples()
    groupingCols <- groupingCols()
    req(gps, input$SV_Types, groupingCols)
    gps$nSvs <- 0
    gps$rate <- 0
    for(svType in input$SV_Types) {
        gps$nSvs <- gps$nSvs + gps[[svType]]
        gps$rate <- gps$rate + gps[[svType]] / gps$coverage
    }
    groups <- gps[, .(
        nProjects = length(unique(project)),
        nSamples = .N,
        coverage = sum(coverage),
        nSvs = sum(nSvs),
        meanSvRate = round(mean(rate), 4),
        sdSvRate = round(sd(rate), 4),
        rates = list(rate),
        projects = list(project)
    ), by = groupingCols]
    groups[, groupLabels := sapply(1:nrow(groups), function(i) {
        paste(sapply(1:length(groupingCols), function(j) {
            paste(groupingCols[j], groups[i][[groupingCols[j]]], sep = " = ")
        }), collapse = ", ")
    })]
    conditionOrder(NULL)
    groupOrder(TRUE)
    groups
})
groupsTable <- bufferedTableServer(
    "groups",
    id,
    input,
    reactive({ 
        groups <- groups() 
        cols <- names(groups)
        groups[, .SD, .SDcols = cols[!(cols %in% c("rates","groupLabels"))]]
    }),
    # editBoxes = list(),
    selection = 'none',
    options = list(
        paging = FALSE,
        searching = FALSE
    )
)

#----------------------------------------------------------------------
# plot element ordering
#----------------------------------------------------------------------
getBucketList <- function(rankListId, labels){
    keepRankListId <- paste0(rankListId, "Keep")
    dropRankListId <- paste0(rankListId, "Drop")
    rankListOnChangeId <- session$ns(paste0(rankListId, "OnChange"))
    bucket_list(
        "",
        add_rank_list(
            text = "Show on Plot",
            labels = labels,
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
            text = "Omit from Plot",
            labels = character(),
            input_id = dropRankListId,
            css_id = dropRankListId,
            options = sortable_options(
                group = rankListId,
                multiDrag = TRUE
            ),
            class = "default-sortable" 
        )
    ) 
}
#----------------------------------------------------------------------
conditionOrder <- reactiveVal(NULL)
output$conditionOrder <- renderUI({
    groupingCols <- groupingCols()
    req(groupingCols)
    getBucketList("conditionOrder", groupingCols)
})
observeEvent(input$conditionOrderOnChange, {
    groupingCols <- groupingCols()
    conditionOrder( sapply(input$conditionOrderOnChange, function(x) which(groupingCols == x)) ) 
})
#----------------------------------------------------------------------
groupOrder <- reactiveVal(TRUE)
output$groupOrder <- renderUI({
    groups <- groups()
    req(groups)
    getBucketList("groupOrder", groups$groupLabels)
})
observeEvent(input$groupOrderOnChange, {
    groups <- groups()
    groupOrder( sapply(input$groupOrderOnChange, function(x) which(groups$groupLabels == x)) ) 
})

#----------------------------------------------------------------------
# output plot
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
ratePlot <- staticPlotBoxServer(
    "ratePlot",
    maxHeight = "400px",
    points  = TRUE, # set to TRUE to expose relevant plot options
    title   = TRUE,
    settings = NULL, # an additional settings template as a list()
    # ... # additional arguments passed to settingsServer
    size = "m",
    create = function() {
        groups <- groups()
        groupingCols <- groupingCols()
        groupOrder <- groupOrder()
        conditionOrder <- conditionOrder()
        if(is.null(conditionOrder)) conditionOrder <- 1:length(groupingCols)
        groups <- groups[groupOrder]
        uniqueProjects <- unique(unlist(groups$projects))
        colors <- CONSTANTS$plotlyColors[1:length(uniqueProjects)]
        names(colors) = uniqueProjects
        maxRate <- max(unlist(groups$rates), na.rm = TRUE)
        nGroups <- nrow(groups)
        nConditions <- length(conditionOrder)
        req(nGroups > 0, nConditions > 0)
        par(mar = c(0.1 + nConditions + 1, 8.1, 0.1, 0.1))
        ratePlot$initializeFrame(
            xlim = c(0.5, nGroups + 0.5),
            ylim = c(0, maxRate) * 1.05,
            xlab = "",
            ylab = "# SVs / Target Coverage",
            yaxs = "i",
            xaxt = "n"
        )
        rect( # make the bar plot
            1:nGroups - 0.35, 
            0, 
            1:nGroups + 0.35, 
            groups[, meanSvRate], 
            lty = 1, 
            lwd = 1,
            col = "grey80"
        )        
        for(i in 1:nGroups){ # overplot individual data points on the bar plot
            rates <- unlist(groups[i, rates])
            projects <- unlist(groups[i, projects])
            ratePlot$addPoints(
                x = jitter2(rates, i - 0.25, i + 0.25),
                y = rates,
                col = sapply(projects, function(x) colors[[x]])
            )            
        }
        mtext( # labels for the condition grid rows
            gsub("_", " ", groupingCols[conditionOrder]), 
            side = 1, 
            line = 1:nConditions, 
            at = 0,
            adj = 1,
                cex = 0.9
        )
        for(i in 1:nConditions){ # fill the condition grid with values by row
            j <- conditionOrder[i]
            mtext(
                unlist(groups[, .SD, .SDcols = groupingCols[j]]), 
                side = 1, 
                line = i, 
                at = 1:nGroups
            )
        }
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
    savedPlots$list  <- bm$outcomes$plots
    savedPlots$names <- bm$outcomes$plotNames
    selfDestruct$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(
        plots     = reactive(savedPlots$list),
        plotNames = reactive(savedPlots$names)
    ),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
