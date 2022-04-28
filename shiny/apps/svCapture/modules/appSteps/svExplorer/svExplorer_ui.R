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

        # selecting sample sources
        fluidRow( box(
            width = 12,
            sampleSelectorUI(ns('sampleSelector'))
        ) ),
        
        # composite views of all filtered SVs
        fluidRow(
            staticPlotBoxUI(
                ns('svLocations'),
                width = 6,
                title = "SV Locations"
            ),
            staticPlotBoxUI(
                ns('svProperties'),
                width = 6,
                title = "SV Properties"
            )
        ),
        
        # sortable table of SVs
        fluidRow(
            box(
                width = 12,
                bufferedTableUI(ns('svsTable'))
            )
        ),
        
        # map image view of a single junction
        fluidRow(
            box(
                width = 12,                
                collapsible = TRUE,
                collapsed = FALSE,
                title = tagList(
                    "Junction Evidence Map",
                    settingsUI(ns('mapSettings'))
                ),
                imageOutput(ns("junctionMapImage"), inline = TRUE)
            )
        ),

        fluidRow(

            # plot of node points of all molecules matching a specific selected junction
            staticPlotBoxUI(
                ns('junctionNodes'),
                width = 4,
                title = "Supporting Molecules"
            ),
            
            # base-level text view of the junction
            box(
                width = 8,
                collapsible = TRUE,
                collapsed = FALSE,
                title = tagList(
                    "Junction Alignment to Reference Genome",
                    settingsUI(ns('alignmentSettings'))
                ),
                tags$style(HTML("
                    .junction { background-color: #ddd;}
                    .alignment .base_A { color: green; }
                    .alignment .base_C { color: blue; }
                    .alignment .base_G { color: brown; }
                    .alignment .base_T { color: red; }
                    .referenceGenome { font-style: oblique; }
                ")),
                div(
                    htmlOutput(ns("junctionAlignment")),
                    style = "font-family: monospace; font-weight: bold; width: 100; overflow: auto;"
                )
            )            
        ),
        ""
    )    
}
