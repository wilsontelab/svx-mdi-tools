#----------------------------------------------------------------------
# UI components for the annotateCells appStep module
#----------------------------------------------------------------------

# module ui function
annotateCellsUI <- function(id, options) {

    # initialize namespace
    module <- 'annotateCells'
    appStepDir <- getAppStepDir(module)
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$annotateCells)

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
        settings = FALSE,

        # track browser styles
        tags$style(slurpFile(file.path(appStepDir, paste0(module, ".css")))),
        # tags$script(slurpFile(file.path(appStepDir, paste0(module, ".js")))),

        # table to select a single scCNV sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 12, 
            collapsible = TRUE
        ),

        # row of cell annotation inputs
        tags$div(
            style = "white-space: nowrap;",
            tags$div(
                class = "annotateCellsInput",
                style = "width: 300px",
                radioButtons(ns("windowSize"), "Window Size (# of 20 kb bins)",
                             choices = 2 ^ (0:6), selected = 1, inline = TRUE)
            ),
            tags$div(
                class = "stepperInput",
                listStepperButtonsUI(ns("cellStepper"))
            ),
        ),
        
        # row of plots
        fluidRow(
            
            # left column of boxes for GC bias plots
            column(width = 3,
                # style = unpad,  
                interactiveScatterplotUI(ns('NR_map_w_vs_gc'), height = '250px'),         
                ""                  
            ),
            
            # right column of boxes for genomic bin plots
            column(width = 9,
                # style = unpad,
                interactiveScatterplotUI(ns('NR_map_w_vs_bin'),  height = '250px')              
            )
        ),

        # appStep UI elements, populate as needed
        ""
    )
}
