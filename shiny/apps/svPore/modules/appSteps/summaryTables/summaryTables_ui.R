#----------------------------------------------------------------------
# UI components for the summaryTables appStep module
#----------------------------------------------------------------------

# module ui function
summaryTablesUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$summaryTables)

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
        svPoreSampleSelectorUI(ns("data")),
        fluidRow(
            bufferedTableBoxUI(
                ns("junctionEdges"),
                title = "Junction Edges",
                width = 6
            ),    
            bufferedTableBoxUI(
                ns("junctionClusters"),
                title = "Junction Clusters",
                width = 6     
            )
        ),
        fluidRow(
            bufferedTableBoxUI(
                ns("segments"),
                title = "Segments",
                width = 6
            ),    
            bufferedTableBoxUI(
                ns("segmentClusters"),
                title = "Segment Clusters",
                width = 6     
            )
        )
    )
}
