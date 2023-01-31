#----------------------------------------------------------------------
# UI components for the collateCNVs appStep module
#----------------------------------------------------------------------

# module ui function
collateCNVsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$collateCNVs)

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
        settings = TRUE,

        # box for selecting data source
        dataSourceTableUI(
            ns("source"), 
            "Project Source", 
            width = 6, 
            collapsible = TRUE
        ),

        fluidRow(
            bufferedTableUI(
                ns("cnvsTable"), 
                title = NULL, 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE
            )            
        ),
        
        fluidRow( box(
            width = 12,
            uiOutput(ns("cellPlot"))            
        )),

        NULL
    )
}
