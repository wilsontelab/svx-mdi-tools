#----------------------------------------------------------------------
# UI components for the normalizeBatch appStep module
#----------------------------------------------------------------------

# module ui function
normalizeBatchUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalizeBatch)

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

        # table to select a single scCNV sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 12, 
            collapsible = TRUE
        ),
        
        # the main action button
        fluidRow(
            column(
                width = 12,
                uiOutput(ns("executeUI"))
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("windowSize"), 
                "Window Size",
                # ...,   
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER
            ),
            staticPlotBoxUI(
                ns("batchEffect"), 
                "Batch Effect",
                # ...,   
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER
            )
        ),
        ""
    )
}
