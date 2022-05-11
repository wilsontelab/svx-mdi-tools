#----------------------------------------------------------------------
# static components to plot summary results over all selected samples
#----------------------------------------------------------------------

# module ui function
analyzeInsertionsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$analyzeInsertions)

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
        
        # composite views of all filtered SVs
        fluidRow(
            staticPlotBoxUI(
                ns('svCounts'),
                width = 6,
                title = "# of SVs By Insert Size"
            ),
            staticPlotBoxUI(
                ns('templateYield'),
                width = 6,
                title = "Yield of Found Templates"
            ),            
            staticPlotBoxUI(
                ns('templateLocations'),
                width = 12,
                title = "Templated Insertion Locations"
            )
        ),

        # sortable table of SVs
        fluidRow(
            filteredSvsTableUI(ns, width = 12)
        )
        # ,

        # # dotplot alignments of junction and genome reference sequences
        # fluidRow(
        #     staticPlotBoxUI(
        #         ns('jxn_side1'),
        #         width = 6,
        #         title = "Junction vs. Left Side"
        #     ),
        #     staticPlotBoxUI(
        #         ns('jxn_side2'),
        #         width = 6,
        #         title = "Junction vs. Right Side"
        #     )
        # )
    )    
}
