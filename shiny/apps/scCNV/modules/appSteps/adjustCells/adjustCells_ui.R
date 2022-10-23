#----------------------------------------------------------------------
# UI components for the adjustCells appStep module
#----------------------------------------------------------------------

# module ui function
adjustCellsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    module <- 'adjustCells'
    appStepDir <- getAppStepDir(module)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$adjustCells)

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
        tags$style(slurpFile(file.path(appStepDir, "adjustCells.css"))),
        tags$script(slurpFile(file.path(appStepDir, "adjustCells.js"))),
        tags$div(
            style = "margin-top: 0; margin-bottom: 10px;",
            tags$div(
                class = "cellPageInput",
                radioButtons(ns('cellType'), "Cell Type", 
                             choices = c("Keep", "Reject"), inline = TRUE,
                             width = "125px")
            ),
            tags$div(
                class = "cellPageInput",
                numericInput(ns('modalCN'), "Modal CN", value = NULL, width = "100px"),
            ),
            tags$div(
                class = "cellPageInput",
                numericInput(ns('cellsPerPage'), "Cells Per Page", value = 5, width = "100px"),
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
                actionButton(ns('nextPage'), ">")
            ),
            tags$div(
                style = "margin-left: 15px; margin-top: 30px; display: inline-block;",
                textOutput(ns("nCellsInType"))
            )
        ),

        # stacked individual cell plots
        uiOutput(ns("cellPlots"))    
    )
}
