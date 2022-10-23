#----------------------------------------------------------------------
# UI components for the adjustFit appStep module
#----------------------------------------------------------------------

# module ui function
adjustFitUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    module <- 'adjustFit'
    appStepDir <- getAppStepDir(module)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$adjustFit)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = FALSE,

        # table to select a single scCNV sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 12, 
            collapsible = TRUE
        ),

        # pagination controls
        tags$style(slurpFile(file.path(appStepDir, "adjustFit.css"))),
        tags$script(slurpFile(file.path(appStepDir, "adjustFit.js"))),
        tags$div(
            style = "margin-bottom: 15px;",
            tags$div(
                class = "cellPageInput",
                selectInput(ns('cellType'), "Cell Type", 
                             choices = c("autoKeep", "autoReject", "userKeep", "userReject"), 
                             width = "125px")
            ),
            tags$div(
                class = "cellPageInput",
                numericInput(ns('cellsPerPage'), "Cells Per Page", value = 10, width = "100px"),
            ),
            tags$div(
                class = "cellPageInput",
                actionButton(ns('prevPage'), "<"),
            ),
            tags$div(
                class = "cellPageInput",
                textInput(ns('pageNumber'), "Page", value = 1, width = "41px"),
            ),
            tags$div(
                class = "cellPageInput",
                actionButton(ns('nextPage'), ">"),
            )
        ),

        # stacked individual cell plots
        uiOutput(ns("cellPlots"))    
        # mdiInteractivePlotUI(ns("cellPlots"))
    )
}
