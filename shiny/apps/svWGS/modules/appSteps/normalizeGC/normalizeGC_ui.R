#----------------------------------------------------------------------
# UI components for the normalizeGC appStep module
#----------------------------------------------------------------------

# module ui function
normalizeGCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalizeGC)

    # UI functions
    plotBox_ <- function(title, column1UI, column2UI, column3UI = NULL){
        box(
            title = title,
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            # collapsible = TRUE,
            column(
                width = 3,
                column1UI
            ),
            column(
                width = 9,
                column2UI
            ),
            column(
                width = 12,
                column3UI
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
        # settings = TRUE,

        # data source selectors
        fluidRow(
            dataSourceTableUI(
                ns("source"), 
                "Data Source", 
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            bufferedTableUI(
                ns("sampleGenome"),
                "Sample with Genome",
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = FALSE
            )            
        ),

        # GC fit and chromosome profile plots
        fluidRow(
            plotBox_(
                "GC Bias Fit",
                tags$div(
                    numericInput(ns("ploidy"),     "Median Ploidy", value = 2, min = 1, max = 4, step = 1),
                    selectInput(ns("expectedSex"), "Expected Sex",  selected = "auto", choices = c("auto","XX","XY"))
                ),
                interactiveScatterplotUI(ns("gcBiasPlot"), height = '400px')
            ),
            plotBox_(
                title = "Chromosome Profiles",
                tags$div(
                    style = "max-height: 400px; overflow: auto;",
                    radioButtons(ns("densityChrom"), "Chromosome", choices = c("genome"), inline = FALSE)
                ),
                interactiveScatterplotUI(ns("chromDensityPlot"), height = '400px')
            )
        ),

        NULL
    )
}
