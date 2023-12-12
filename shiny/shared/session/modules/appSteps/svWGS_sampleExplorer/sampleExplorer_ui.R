#----------------------------------------------------------------------
# UI components for the sampleExplorer appStep module
#----------------------------------------------------------------------

# module ui function
svWGS_sampleExplorerUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$svWGS_sampleExplorer)

    # UI functions
    plotBox_ <- function(type, title){
        staticPlotBoxUI(
            ns(paste0(type, "Plot")), 
            title = title,
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            collapsed = FALSE         
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

        # data source selectors
        svWGS_explorer_dataSelectorsUI(ns("dataSelectors")),
        
        # microhomology plots and tables
        fluidRow(
            # box(
            #     title = "Group By",
            #     width = 2,
            #     solidHeader = TRUE,
            #     status = "primary",
            #     collapsible = TRUE,
            #     collapsed = FALSE,
            #     tags$div(
            #         style = "padding-left: 10px; padding-top: 5px;",
            #         checkboxGroupInput(
            #             ns("microhomologyGroupBy"),
            #             NULL,
            #             choices = c("Project", "Sample_ID", "Genome", "Clonal", "SV_Type"),
            #             inline = FALSE
            #         )
            #     )
            # ),
            plotBox_("countCorrelation", "Junction Count Correlation")
        ),
        # fluidRow(
        #     column(width = 2),
        #     plotBox_("sizeCorrelation", "Insert Size vs. SV Size")
        # ),
        NULL
    )
}
