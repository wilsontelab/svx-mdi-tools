#----------------------------------------------------------------------
# UI components for the summaryPlots appStep module
#----------------------------------------------------------------------

# module ui function
summaryPlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$summaryPlots)

    axisColumns <- c("eventSize", "insertSize","nJunctionInstances","hasAdapter","blastIdentity","gapCompressedIdentity")
    colorColumns <- c("edgeType","passedFlankCheck","passedBandwidth","hasAdapter")
    filterTypes <- c("AllJunctions","RealJunctions","MatchableJunctions","FusableJunctions","SingletonJunctions")
    filterLabels <- c("passedFlankCheck","passedBandwidth","hasAdapter","isChrM","invert")
    filterLabelsSelected <- c("passedFlankCheck","passedBandwidth","hasAdapter","isChrM")

# mapQ
# nStrands

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        # settings = TRUE,

        # appStep UI elements, populate as needed
        fluidRow(
            dataSourceTableUI(ns("source"), "Sample", width = 8, collapsible = FALSE)
        ),
        fluidRow(
            # column(width = 2, selectInput(ns("xAxisColumn"), "X Axis", choices = axisColumns, selected = "nJunctionInstances")),
            # column(width = 2, selectInput(ns("yAxisColumn"), "Y Axis", choices = axisColumns, selected = "eventSize")),
            column(width = 2, selectInput(ns("xAxisColumn"), "X Axis", choices = axisColumns, selected = "eventSize")),
            column(width = 2, selectInput(ns("yAxisColumn"), "Y Axis", choices = axisColumns, selected = "insertSize")),
            column(width = 2, selectInput(ns("colorColumn"), "Color By", choices = colorColumns, selected = "edgeType")),
            column(width = 2, numericInput(ns("opacity"), "Point Opacity", 0.5, min = 0.5, max = 1, step = 0.025)),
            column(width = 2, selectInput(ns("edgeType"), "Edge Type", choices = junctionTypes, selected = "All"))
            # column(width = 2, numericInput(ns("minMapQ"), "Min MapQ", 55, min = 0, max = 60, step = 5)),
        ),
        fluidRow(
            style = "margin-bottom: 20px;",
            column(width = 2, radioButtons(ns("showChrM"), "Show ChrM", c("No", "Yes"), inline = TRUE, selected = "No")),             
            column(width = 10, radioButtons(ns("filterType"), "Filter Type", filterTypes, inline = TRUE, selected = "AllJunctions"))
            # column(width = 12, checkboxGroupInput(ns("filters"), "Filters", filterLabels, inline = TRUE, selected = filterLabelsSelected)) 
        ),
        # fluidRow(
        #     style = "margin-bottom: 20px;",
        #     column(width = 12, checkboxGroupInput(ns("filters"), "Filters", filterLabels, inline = TRUE, selected = filterLabelsSelected)) 
        # ),
        fluidRow(
            box(
                width = 12,
                title = "Summary Plot",
                mdiInteractivePlotUI(ns("summaryPlot"))
            )
        ),
        fluidRow(
            box(
                width = 12,
                title = "Molecule Plot",
                # plotOutput(ns("moleculePlot"), height = "600px")
                mdiInteractivePlotUI(ns("moleculePlot"))
            )
        )
        # ,
        # fluidRow(
        #     bufferedTableUI (
        #         ns("junctionsTypesTable"), 
        #         title = "Junction Types", 
        #         downloadable = TRUE,
        #         width = 4,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     ),
        #     staticPlotBoxUI(
        #         ns("svTrianglePlot"), 
        #         title = "Triangle Plot",
        #         width = 4
        #     )
        # ), 
        # fluidRow(
        #     staticPlotBoxUI(
        #         ns("positionDensityPlot"), 
        #         title = "Position Density",
        #         width = 4
        #     ),            
        #     staticPlotBoxUI(
        #         ns("sizeDensityPlot"), 
        #         title = "Size Density",
        #         width = 4
        #     )
        # ),   
        # fluidRow(
        #     bufferedTableUI (
        #         ns("junctionsTable"), 
        #         title = "Junctions", 
        #         downloadable = TRUE,
        #         width = 8,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     ) 
        # )
    )
}
