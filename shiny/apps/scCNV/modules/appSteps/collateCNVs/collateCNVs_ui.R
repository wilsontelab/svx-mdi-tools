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
        collapsibleDivUI(
            ns("leaderText"), 
            "Select one or more project sources that might share the same kept CNVs. Scanning through the chromosomes,",
            "click to select CNVs that appear to be identical. Click 'COMMIT GROUP' to associate them with each other.",
            "Occasionally, you may also group two adjacent CNVs in the same cell, e.g., when you think ",
            "they were inappropriately split by the pipeline.",
            "In all cases, the CNVs being grouped should have the same copy number signature, e.g., representing a single-copy deletion.",
            "CNVs not placed into a group are considered to form their own single-member group."                                
        ),

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
            width = 8, 
            collapsible = TRUE
        ),
 
        fluidRow(
            bufferedTableUI(
                ns("cnvsTable"), 
                title = "CNVs Table", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = TRUE
            )
        ),

        # stacked individual cell plots
        tags$div(
            style = "white-space: nowrap;",
            tags$div(
                id = ns("chromPlotsWrapper"),
                style = "min-width: 995px;",
                class = "cellPlotsWrapper",
                uiOutput(ns("chromPlots")),
                tags$div(class = "cellStackVertical")
            )            
        ),
        tags$div(actionLink(ns("clearChromAssignments"),  "clear all assignments on this chromosome"),   style = "display: inline-block; margin-top: 10px;"),
        tags$div(actionLink(ns("clearSourceAssignments"), "clear all assignments for selected sources"), style = "display: inline-block; margin-left: 20px;"),
        tags$div(actionLink(ns("clearAllAssignments"),    "clear all assignments across all sources"),   style = "display: inline-block; margin-left: 20px;")
    )
}
