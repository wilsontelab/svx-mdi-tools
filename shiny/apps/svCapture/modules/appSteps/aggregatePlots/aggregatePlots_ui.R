#----------------------------------------------------------------------
# static components to plot summary results over all selected samples
#----------------------------------------------------------------------

# module ui function
aggregatePlotsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$aggregatePlots)

    # incorporate options text into templates
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )
    
    # unpad columns with boxes
    unpad <- "padding-left: 0; padding-right: 0;"

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
        
        # composite views of all filtered SVs
        fluidRow(
            staticPlotBoxUI(
                ns('svRates'),
                width = 6,
                title = "SV Rates"
            )
        ),
        
        # table of plotted values
        fluidRow(
            column(
                width = 12,
                style = unpad,
                box(
                    width = 12,
                    bufferedTableUI(ns('ratesTable'))
                )
            )
        )
    )    
}
