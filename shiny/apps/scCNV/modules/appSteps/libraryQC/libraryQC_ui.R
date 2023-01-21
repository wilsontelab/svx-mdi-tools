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

        # box for selecting sample source
        dataSourceTableUI(
            ns("source"), 
            "Sample", 
            width = 12, 
            collapsible = TRUE
        ),

        fluidRow( box(
            width = 12,
            column(
                width = 12,
                textOutput(ns('nFailedLibraries'), inline = TRUE)           
            )
        ) ),

        # box with plots to explore QC metric relationships and click to reject cells
        fluidRow(
            # mdiInteractivePlotUI(ns("nReads_vs_alignRate")),
            staticPlotBoxUI(
                ns("nReads_vs_alignRate"),
                NULL, 
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER,
                terminal = FALSE 
            ),
            staticPlotBoxUI(
                ns("nReads_vs_dupRate"),
                NULL, 
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER,
                terminal = FALSE 
            ),
            staticPlotBoxUI(
                ns("nReads_vs_windowPower"),
                NULL, 
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER,
                terminal = FALSE 
            ),
            staticPlotBoxUI(
                ns("nReads_vs_cnsd"),
                NULL, 
                documentation = serverEnv$IS_DEVELOPER,
                code = serverEnv$IS_DEVELOPER,
                console = serverEnv$IS_DEVELOPER,
                terminal = FALSE 
            )
        )
    )    
}
