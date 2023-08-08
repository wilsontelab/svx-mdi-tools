#----------------------------------------------------------------------
# UI components for the dataExplorer appStep module
#----------------------------------------------------------------------

# module ui function
dataExplorerUI <- function(id, options) {
    ns <- NS(id)
    options <- setDefaultOptions(options, stepModuleInfo$dataExplorer)
    standardSequentialTabItem(
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
        dataSourceTableUI(
            ns("dataSource"), 
            "Data Source", 
            width = 12, 
            collapsible = TRUE
        ),
        fluidRow(
            box(
                width = 6,
                title = "R code to return a data.table",
                status = 'primary',
                solidHeader = TRUE, 
                collapsible = TRUE,
                style = "padding: 0;",
                tags$div(
                    style = "padding: 5px; text-align: center;",
                    bsButton(ns("updateTableData"), "Update Outputs", style = "success")
                ),
                tags$div(
                    id = ns("tableDataEditor"),
                    style = "height: 400px; border-top: 1px solid grey;"
                )                  
            ), 
            box(
                width = 6,
                title = "Source Object Structure",
                status = "primary",
                solidHeader = TRUE, 
                collapsible = TRUE,
                tags$div(
                    style = "padding: 7px; text-align: center;",
                    radioButtons(ns("sourceObjectChoices"), NULL, inline = TRUE, choices = names(options$sourceObjects))
                ),
                tags$div(
                    style = "height: 400px; border-top: 1px solid grey;", 
                    verbatimTextOutput(ns("sourceObjectStructure"))
                )     
            )
        ),
        fluidRow(
            bufferedTableUI(
                ns("tableContents"), 
                title = "Contents of the data.table", 
                downloadable = TRUE, 
                width = 12, 
                status = 'primary',
                solidHeader = TRUE, 
                collapsible = TRUE
            )
        ),
        fluidRow(
            box(
                width = 6,
                title = "R code to return a plot",
                status = 'primary',
                solidHeader = TRUE, 
                collapsible = TRUE,
                style = "padding: 0;",
                tags$div(
                    style = "padding: 5px; text-align: center;",
                    bsButton(ns("updatePlot"), "Update Plot", style = "success")
                ),
                tags$div(
                    id = ns("plotEditor"),
                    style = "height: 400px; border-top: 1px solid grey;"
                )                  
            ), 
            box(
                width = 6,
                title = "Plot",
                status = "primary",
                solidHeader = TRUE,
                style = "padding: 0;",
                tags$div(
                    style = "height: 435px;", 
                    mdiInteractivePlotUI (ns("plot"))
                )     
            )
        )
    )
}
