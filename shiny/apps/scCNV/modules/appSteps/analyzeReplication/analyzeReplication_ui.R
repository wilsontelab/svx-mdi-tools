#----------------------------------------------------------------------
# UI components for the analyzeReplication appStep module
#----------------------------------------------------------------------

# module ui function
analyzeReplicationUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$analyzeReplication)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        collapsibleDivUI(
            ns("leaderText"), 
            "Select one or more project sources expected to represent cells with a common replication timing."                                
        ),

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = TRUE,

        # box for selecting sample source
        scCnvDataSourceTableUI(ns),

        # stacked individual cell plots
        tags$div(
            style = "white-space: nowrap;",
            tags$div(
                id = ns("chromPlotsWrapper"),
                class = "cellPlotsWrapper",
                uiOutput(ns("chromPlots")),
                tags$div(class = "cellStackVertical")
            )            
        ),
        NULL    
    )
}
