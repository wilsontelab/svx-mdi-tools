#----------------------------------------------------------------------
# UI components for the sampleSummary appStep module
#----------------------------------------------------------------------

# module ui function
sampleSummaryUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$sampleSummary)

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
        dataSourceTableUI(
            ns("dataSource"), 
            "Data Source", 
            width = 12, 
            collapsible = TRUE
        ),
        fluidRow(
            box(
                width = 12,
                title = "Junction Networks",
                status = 'primary',
                solidHeader = TRUE, 
                collapsible = FALSE,
                style = "padding: 0; height: 800px;",
                imageOutput(ns("graphImage"))              
            )
        )
    )
}
