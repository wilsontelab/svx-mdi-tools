#----------------------------------------------------------------------
# UI components for the keepRejectCells appStep module
#----------------------------------------------------------------------

# module ui function
keepRejectCellsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    module <- 'keepRejectCells'
    appStepDir <- getAppStepDir(module)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$keepRejectCells)

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

        # box for selecting data source
        dataSourceTableUI(
            ns("source"), 
            "Project Source", 
            width = 6, 
            collapsible = TRUE
        ),

        # pagination controls
        includeCSS(file.path(appStepDir, "keepRejectCells.css")),
        tags$head(includeScript(file.path(appStepDir, "keepRejectCells.js"))),
        tags$div(
            style = "margin-top: 0; margin-bottom: 10px; white-space: nowrap;",
            tags$div(
                class = "cellPageInput",
                style = "width: 175px; margin-right: 10px;", 
                selectInput(ns('sampleNameFilter'), "Sample", choices = c(All = "__ALL__"), width = "100%")
            ),
            tags$div(
                class = "cellPageInput",
                style = "width: 175px; margin-right: 10px;", 
                selectInput(ns("cellIdFilter"), "Cell", choices = c(All = "__ALL__"), width = "100%")
            ),
            tags$div(
                class = "cellPageInput",
                radioButtons(ns('cellStatus'), "Cell Status", 
                             choices = c("All", "Keep", "Reject", "Bad"), inline = TRUE,
                             width = "235px")
            ),
            tags$div(
                class = "cellPageInput",
                radioButtons(ns('replicating'), "Replicating", 
                             choices = c("All", "Yes", "No"), inline = TRUE,
                             width = "150px")
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
                style = "margin-left: 15px; margin-top: 30px; display: inline-block; white-space: nowrap;",
                textOutput(ns("nCellsInType"),inline = TRUE),
                actionLink(ns("clearImageCache"), "clear image cache", style = "margin-left: 10px;"),
                actionLink(ns("resetEverything"), "reset everything", style = "margin-left: 10px;")
            )
        ),

        # stacked individual cell plots
        tags$div(
            style = "white-space: nowrap;",
            tags$div(
                id = ns("genomePlotsWrapper"),
                class = "cellPlotsWrapper",
                # uiOutput(ns("genomeLabelRow")),
                uiOutput(ns("genomePlots")),
                tags$div(class = "cellStackVertical")
            ),
            tags$div(
                id = ns("chromPlotsWrapper"),
                class = "cellPlotsWrapper",
                uiOutput(ns("chromPlots")),
                tags$div(class = "cellStackVertical")
            )            
        )
    )
}
