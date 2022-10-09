#----------------------------------------------------------------------
# UI components for the markCells appStep module
#----------------------------------------------------------------------

# module ui function
markCellsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    module <- 'markCells'
    appStepDir <- getAppStepDir(module)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$markCells)

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
        settings = TRUE,

        # table to select a single scCNV sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 12, 
            collapsible = TRUE
        ),

        # pagination controls
        tags$style(slurpFile(file.path(appStepDir, "markCells.css"))),
        tags$div(
            style = "margin-bottom: 15px;",
            tags$div(
                class = "cellPageInput",
                radioButtons(ns('Cell_Type'), "", choices = c("good", "bad"), inline = TRUE, width = "125px")
            ),
            tags$div(
                class = "cellPageInput",
                numericInput(ns('Cells_Per_Page'), "Cells Per Page", value = 10, width = "100px"),
            ),
            tags$div(
                class = "cellPageInput",
                actionButton(ns('prevPage'), "<"),
            ),
            tags$div(
                class = "cellPageInput",
                textInput(ns('Page_Number'), "Page", value = 1, width = "41px"),
            ),
            tags$div(
                class = "cellPageInput",
                actionButton(ns('nextPage'), ">"),
            )
        ),

        # stacked individual cell plots        
        mdiInteractivePlotUI(ns("cellPlots"))
    )
}
