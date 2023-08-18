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
            mdiInteractivePlotUI(ns("pairedQualityPlot"))
        ),
        fluidRow(
            style = "padding: 5px;",
            column(
                offset = 2,
                width = 2,
                style = "padding-top: 7px;",
                radioButtons(ns("moleculeTypeFilter"), NULL, choices = c("Kept", "Rejected"), selected = "Kept", inline = TRUE)
            ),
            column(
                width = 4,
                listStepperButtonsUI(ns("moleculeTypeStepper")) 
            )
        ),             
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
        )
    )
}
