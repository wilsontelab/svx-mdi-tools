#----------------------------------------------------------------------
# static components to filter and examine SV locations and junction sequences
#----------------------------------------------------------------------

# module ui function
svExplorerUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$svExplorer)

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
        
        # # composite views of all filtered SVs
        # fluidRow(
        #     staticPlotBoxUI(
        #         ns('svLocations'),
        #         width = 6,
        #         title = "SV Locations"
        #     ),
        #     staticPlotBoxUI(
        #         ns('svProperties'),
        #         width = 6,
        #         title = "SV Properties"
        #     )
        # ),
        
        # box with sortable table of SVs
        fluidRow(
            box(
                width = 12,
                bufferedTableUI(ns('svsTable'))
            )
        ),

        fluidRow(

            # # box with the node points of all molecules matching a specific selected junction
            # staticPlotBoxUI(
            #     ns('junctionNodes'),
            #     width = 4,
            #     title = "Supporting Molecules"
            # ),

            # # box with image view of a single junction
            # box(
            #     width = 8,
            #     imageOutput(ns("junctionMapImage")),
            # )

            # box(
            #     width = 8,
            #     span(
            #         htmlOutput(ns("junctionAlignment")),
            #         style = "font-family: monospace; font-weight: bold;"
            #     )
            # )            
        ),
        
        # box with text view of a single junction
        fluidRow(
            box(
                width = 12,
                span(
                    htmlOutput(ns("junctionAlignment")),
                    style = "font-family: monospace; font-weight: bold;"
                )
            )    
        ),
        ""
    )    
}
