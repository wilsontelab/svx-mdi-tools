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

        # page header links
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        # settings = TRUE,

        # table to select a single scCNV sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 12, 
            collapsible = TRUE
        ),

        # summary stats and plots
        fluidRow(
            box(
                width = 12,
                title = "Cell Counts",
                uiOutput(ns("metrics"))
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("windowSizes"), 
                "Window Sizes",
                # ...,   
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER
            )
        )
    )
}
