#----------------------------------------------------------------------
# UI components for the assemblyPlots appStep module
#----------------------------------------------------------------------

# module ui function
assemblyPlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$assemblyPlots)

    # bucket and static plot box
    assemblyPlotBox <- function(type, title, collapsed = TRUE){
        fluidRow(
            box(
                width = 12,
                title = title,
                collapsible = TRUE,
                collapsed = collapsed,
                fluidRow(
                    column(
                        width = 6,
                        style = "padding: 0;",
                        uiOutput(ns(paste("conditions", type, sep = "_"))),
                        uiOutput(ns(paste("groups", type, sep = "_")))
                    ),
                    staticPlotBoxUI(
                        ns(paste0(type, "Plot")), 
                        title = "",
                        width = 6,
                        collapsible = TRUE,
                        collapsed = collapsed
                    )
                )
            )
        )
    }

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
                ns("savedPlotSets"), 
                title = "Saved Plot Assemblies", 
                width = 10,
                collapsible = TRUE,
                collapsed = FALSE,
                downloadable = TRUE 
            ),
            column(
                width = 2,
                bsButton(
                    ns("savePlotSet"),
                    "Save Plots",
                    style = "success",
                    block = TRUE
                )
            ),
            column(
                width = 12,
                uiOutput(ns('savePlotSetFeedback'))
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
                    tags$div(
                        width = 3,
                        bsButton(
                            ns("suspendDataProcessing"),
                            "Allow Data Processing",
                            style = "warning",
                            block = TRUE,
                            type = "toggle",
                            value = TRUE
                        )                       
                    ),
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
                collapsed = TRUE,
                downloadable = TRUE
            )
        ),

        # output plots
        assemblyPlotBox("svFrequencies", "SV Frequencies", collapsed = FALSE),
        assemblyPlotBox("microhomology", "Microhomology/Insert Distributions", collapsed = TRUE),
        assemblyPlotBox("endpoints", "SV Endpoint Distributions", collapsed = TRUE),
        assemblyPlotBox("svSizes", "SV Size Distributions", collapsed = TRUE),
        # assemblyPlotBox("insertionTemplates", "Insertion Templates", collapsed = TRUE),
        # assemblyPlotBox("junctionBases", "Junction Base Usage", collapsed = TRUE),
        NULL
    )
}
