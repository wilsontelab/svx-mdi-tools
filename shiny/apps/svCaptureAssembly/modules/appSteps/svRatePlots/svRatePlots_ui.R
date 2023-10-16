#----------------------------------------------------------------------
# UI components for the svRatePlots appStep module
#----------------------------------------------------------------------

# module ui function
svRatePlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$svRatePlots)

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

        # saving figures
        fluidRow(
            bufferedTableUI(
                ns("savedPlots"), 
                title = "Saved Plots", 
                width = 10,
                collapsible = TRUE,
                collapsed = FALSE,
                downloadable = TRUE 
            ),
            column(
                width = 2,
                bsButton(
                    ns("savePlot"),
                    "Save Plot",
                    style = "success",
                    block = TRUE
                )
            ),
            column(
                width = 12,
                uiOutput(ns('savePlotFeedback'))
            )
        ),

        # assembly selection and display of all its samples
        dataSourceTableUI(
            ns("dataSource"), 
            "Assembly Set", 
            width = 12, 
            collapsible = TRUE
        ),
        fluidRow(
            bufferedTableUI(
                ns("allSamples"),
                title = "All Assembled Samples", 
                width = 12,
                collapsible = TRUE,
                collapsed = TRUE,
                downloadable = TRUE
            )
        ),

        # selection of columns and values to plot
        fluidRow(
            box(
                title = "Plot Settings", 
                width = 12,
                tags$div(
                    style = "padding: 10px;",
                    checkboxGroupInput(
                        ns("Group_By"),
                        "Group By",
                        choices = NULL,
                        selected = NULL,
                        inline = TRUE,
                        width = "100%"
                    ),
                    checkboxGroupInput(
                        ns("Required"),
                        "Required",
                        choices = NULL,
                        selected = NULL,
                        inline = TRUE,
                        width = "100%"
                    ),
                    checkboxGroupInput(
                        ns("Prohibited"),
                        "Prohibited",
                        choices = NULL,
                        selected = NULL,
                        inline = TRUE,
                        width = "100%"
                    ),
                    checkboxGroupInput(
                        ns("SV_Types"),
                        "SV Types",
                        choices = c("deletion", "duplication", "inversion", "translocation"),
                        selected = "deletion",
                        inline = TRUE,
                        width = "100%"
                    ),
                    checkboxGroupInput(
                        ns("Projects"),
                        "Projects",
                        choices = character(),
                        inline = TRUE,
                        width = "100%"
                    )
                )
            )
        ), 

        # tabular views of the matching samples and groups
        fluidRow(
            bufferedTableUI(
                ns("groupedProjectSamples"), 
                title = "Matching Samples", 
                width = 12,
                collapsible = TRUE,
                collapsed = TRUE,
                downloadable = TRUE
            )
        ),
        fluidRow(
            bufferedTableUI(
                ns("groups"), 
                title = "Groups to Plot", 
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE,
                downloadable = TRUE
            )
        ),

        # output plot and item ordering
        fluidRow(
            column(
                width = 6,
                style = "padding: 0;",
                uiOutput(ns("conditionOrder")),
                uiOutput(ns("groupOrder"))                
            ),
            staticPlotBoxUI(
                ns("ratePlot"), 
                "SV Rate Plot",
                width = 6 # additional arguments passed to shinydashboard::box()    
            )
        )
    )
}
