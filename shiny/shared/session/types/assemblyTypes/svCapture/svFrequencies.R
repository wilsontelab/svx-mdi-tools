
assemblyPlotBoxServer <- function(id, session, input, output, 
                                  isProcessingData, groupingCols, groups){
    condId <- paste("conditions", id, sep = "_")
    condChangeId <- paste0(condId, "OnChange")
    conditionsReactive <- reactiveVal()
    output[[condId]] <- renderAssemblyConditionsBucket(session, groupingCols, id, conditionsReactive)
    observeEvent(input[[condChangeId]], { conditionsReactive(input[[condChangeId]]) })
    #----------------------------------------------------------------------
    grpId <- paste("groups", id, sep = "_")
    grpChangeId <- paste0(grpId, "OnChange")
    groupsReactive <- reactiveVal()
    output[[grpId]] <- renderAssemblyGroupsBucket(session, groups, id, groupsReactive)
    observeEvent(input[[grpChangeId]], { groupsReactive(input[[grpChangeId]]) })
    #----------------------------------------------------------------------
    data <- reactive({
        req(isProcessingData())
        conditions_ <- conditionsReactive()
        groups_ <- groupsReactive()
        req(conditions_, groups_)
        nConditions <- length(conditions_)
        nGroups     <- length(groups_)

        mar <- APC[[id]]$mar
        mar[1] <- mar[1] + nConditions

    })
    plot <- plotBoxServer
}

svCapture_svFrequencies_plotFrame <- function(){
    list(
        Plot_Frame = getAssemblyPlotFrame(
            plot, 
            plot$settings$get("Groups","Group_Width_Inches") * nGroups, 
            APC[[id]]$insideHeight, 
            mar
        ),
        mar = mar,
        conditions_ = conditions_,
        groups_ = groups_,
        nConditions = nConditions,
        nGroups = nGroups
    )
}

#----------------------------------------------------------------------
# svFrequencies plot
#----------------------------------------------------------------------
# svFrequenciesType <- "svFrequencies"
# conditions_svFrequencies <- reactiveVal()
# output$conditions_svFrequencies <- renderAssemblyConditionsBucket(session, groupingCols, svFrequenciesType, conditions_svFrequencies)
# observeEvent(input$conditions_svFrequenciesOnChange, { conditions_svFrequencies(input$conditions_svFrequenciesOnChange) })
# #----------------------------------------------------------------------
# groups_svFrequencies <- reactiveVal(NULL)
# output$groups_svFrequencies <- renderAssemblyGroupsBucket(session, groups, svFrequenciesType, groups_svFrequencies)
# observeEvent(input$groups_svFrequenciesOnChange, { groups_svFrequencies(input$groups_svFrequenciesOnChange) })
#----------------------------------------------------------------------
svFrequenciesPlotTrigger <- reactiveVal(1)
svFrequenciesData <- reactive({ # not a slow step, not worth caching
    # req(isProcessingData())
    # conditions_ <- conditions_svFrequencies()
    # groups_ <- groups_svFrequencies()
    # req(conditions_, groups_)
    # nConditions <- length(conditions_)
    # nGroups     <- length(groups_)
    # mar <- APC$svFrequencies$mar
    # mar[1] <- mar[1] + nConditions
    # list(
    #     Plot_Frame = getAssemblyPlotFrame(
    #         svFrequenciesPlot, 
    #         svFrequenciesPlot$settings$get("Groups","Group_Width_Inches") * nGroups, 
    #         APC$svFrequencies$insideHeight, 
    #         mar
    #     ),
    #     mar = mar,
    #     conditions_ = conditions_,
    #     groups_ = groups_,
    #     nConditions = nConditions,
    #     nGroups = nGroups
    # )
})
svFrequenciesPlotFrame <- plotFrameReactive(svFrequenciesData) 
svFrequenciesPlot <- staticPlotBoxServer(
    "svFrequenciesPlot",
    settings = assemblyPlotSettings$svFrequencies, 
    size = "m",
    Plot_Frame = svFrequenciesPlotFrame,
    create = function() {
        req(isProcessingData())
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
