#----------------------------------------------------------------------
# static components to show plots with QC values and distributions
# for all samples in a sampleSet, to explore run quality and reject cells
#----------------------------------------------------------------------

# module ui function
libraryQCUI <- function(id, options) {
    ns <- NS(id)    
    options <- setDefaultOptions(options, stepModuleInfo$libraryQC)
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )
    
    # return the UI contents
    standardSequentialTabItem(
        # HTML(paste( options$longLabel, settingsUI(ns('settings')) )),
        HTML(options$longLabel),
        leaderText,

        # page header links
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = TRUE,

        # box for selecting sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 6, 
            collapsible = TRUE
        ),

        # high-level count summaries
        fluidRow( box(
            width = 12,
            uiOutput(ns("sourceSummary"))
        ) ),

        # box with plots to explore QC metric relationships and click to reject cells
        fluidRow( box(
            width = 12,
            mdiInteractivePlotUI(ns("nReads_vs_alignRate")),
            mdiInteractivePlotUI(ns("nReads_vs_dupRate")),
            mdiInteractivePlotUI(ns("nReads_vs_windowPower")),
            mdiInteractivePlotUI(ns("nReads_vs_cnsd")),
        )),
        fluidRow( box(
            width = 12,
            uiOutput(ns("cellPlot"))            
        ))
    )    
}
