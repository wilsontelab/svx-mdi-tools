#----------------------------------------------------------------------
# UI components for the keepReject appStep module
#----------------------------------------------------------------------

# module ui function
keepRejectUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    module <- 'keepReject'
    appStepDir <- getAppStepDir(module)

    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$keepReject)

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
        settings = FALSE,

        # includeCSS(file.path(appStepDir, "keepReject.css")),
        # tags$head(includeScript(file.path(appStepDir, "keepReject.js"))),

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
        # fluidRow(
        #     bufferedTableUI (
        #         ns("junctionsTable"), 
        #         title = "Junctions", 
        #         downloadable = TRUE,
        #         width = 12,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     )
        # ),        
        # fluidRow(
        #     bufferedTableUI (
        #         ns("moleculeTypesTable"), 
        #         title = "Molecule Types", 
        #         downloadable = TRUE,
        #         width = 12,
        #         collapsible = TRUE,
        #         collapsed = FALSE
        #     )
        # ),        
        # fluidRow(
        #     staticPlotBoxUI(
        #         ns("svTrianglePlot"), 
        #         title = "Triangle Plot",
        #         width = 4
        #     ),            
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
        fluidRow(
            box(
                width = 12,
                tags$div(
                    style = "display: inline-block; vertical-align: top;",
                    uiOutput(ns("moleculeMetadata"))
                ),
                tags$div(
                    style = "display: inline-block; width: 800px; vertical-align: top;",
                    plotOutput(ns("moleculeQcPlot"), width = "100%", height ="600px")
                )
            )
            # staticPlotBoxUI(
            #     ns("moleculeDotPlot"), 
            #     title = "Alignment Dot Plot",
            #     width = 12
            # )
        ),
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
