#----------------------------------------------------------------------
# static components to plot summary results for non-SV artifacts
#----------------------------------------------------------------------

# module ui function
artifactExplorerUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$artifactExplorer)

    # incorporate options text into templates
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )

    # return the UI contents
    standardSequentialTabItem(
        HTML(paste( options$longLabel, settingsUI(ns('settings')) )),
        # HTML(options$longLabel),
        leaderText,

        # box for selecting sample sources
        fluidRow( box(
            width = 12,
            sampleSelectorUI(ns('sampleSelector'))
        ) ),
        
        # figures related to t- SVs
        fluidRow(
            staticPlotBoxUI(
                ns('distancePlot'),
                width = 6,
                title = "SV 2nd End Distance"
            ),
            staticPlotBoxUI(
                ns('microhomologyPlot'),
                width = 6,
                title = "Artifact Microhomology"
            ) 
        )
    )    
}
