#----------------------------------------------------------------------
# static components to generate stacked genomic plots and heatmaps
#----------------------------------------------------------------------

# module ui function
compositePlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$compositePlots)

    # incorporate options text into templates
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )
    
    # return the UI contents
    standardSequentialTabItem(
        HTML(paste( options$longLabel, settingsUI(ns('settings')) )),
        leaderText,
        
        # top box for selecting and labeling sample sources
        fluidRow( box(
            
            # select the Sample Set, Group and Type (i.e., set of related samples/cells)
            width = 12,
            column(
                width = 12,
                sampleSetGroupTypeUI(ns('data'))
            ),
            
        )),
        
        # bottom box with genomic plots
        fluidRow(
            box(
                width = 12,
                uiOutput(ns('beCorrectedDepthPlotUI'))
            ),         
        )   
    )    
}
