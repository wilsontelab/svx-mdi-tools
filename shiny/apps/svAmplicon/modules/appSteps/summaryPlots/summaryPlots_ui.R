#----------------------------------------------------------------------
# UI components for the summaryPlots appStep module
#----------------------------------------------------------------------

# module ui function
summaryPlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$summaryPlots)

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
            box(
                width = 4,
                sampleSetUI (ns("sampleSet"))
            ),
            bufferedTableUI (
                ns("ampliconsTable"), 
                title = "Amplicons", 
                downloadable = TRUE,
                width = 8,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),
        fluidRow(
            bufferedTableUI (
                ns("pathClassesTable"), 
                title = "Molecule Path Classes", 
                downloadable = TRUE,
                width = 4,
                collapsible = TRUE,
                collapsed = FALSE
            ),
            bufferedTableUI (
                ns("junctionsTypesTable"), 
                title = "Junction Types", 
                downloadable = TRUE,
                width = 4,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("svCirclesPlot"), 
                title = "Circles Plot",
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ), 
        fluidRow(
            # staticPlotBoxUI(
            #     ns("svTrianglePlot"), 
            #     title = "Triangle Plot",
            #     width = 4,
            #     collapsible = TRUE,
            #     collapsed = TRUE
            # ),        
            staticPlotBoxUI(
                ns("positionDensityPlot"), 
                title = "Position Density",
                width = 6,
                collapsible = TRUE,
                collapsed = TRUE
            ),            
            staticPlotBoxUI(
                ns("sizeDensityPlot"), 
                title = "Size Density",
                width = 6,
                collapsible = TRUE,
                collapsed = TRUE
            )
        ),   
        fluidRow(
            bufferedTableUI (
                ns("junctionsTable"), 
                title = "Junctions", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            ) 
        )
        #     bufferedTableUI (
        #         ns("moleculeTypesTable"), 
        #         title = "Molecule Types", 
        #         downloadable = TRUE,
        #         width = 12,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     )   
        # fluidRow(
        #     bufferedTableUI (
        #         ns("moleculeTypeExpansion"), 
        #         title = "Molecule Type Details", 
        #         downloadable = TRUE,
        #         width = 12,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     )
        # )
    )
}
