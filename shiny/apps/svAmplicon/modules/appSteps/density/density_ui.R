#----------------------------------------------------------------------
# UI components for the density appStep module
#----------------------------------------------------------------------

# module ui function
densityUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$density)

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
        settings = TRUE,

        # appStep UI elements, populate as needed
        fluidRow(box(
            width = 4,
            sampleSetUI (ns("sampleSet"))
        )),
        fluidRow(
            bufferedTableUI (
                ns("ampliconsTable"), 
                title = "Amplicons", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
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
        ),        
        fluidRow(
            bufferedTableUI (
                ns("moleculeTypesTable"), 
                title = "Molecule Types", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),        
        fluidRow(
            staticPlotBoxUI(
                ns("positionDensityPlot"), 
                title = "Position Density"
            ),            
            staticPlotBoxUI(
                ns("svTrianglePlot"), 
                title = "Triangle Plot"
            ),
            staticPlotBoxUI(
                ns("sizeDensityPlot"), 
                title = "Size Density"
            )
        ),        
        fluidRow(
            bufferedTableUI (
                ns("moleculeTypeExpansion"), 
                title = "Molecule Type Details", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )
        )
    )
}
