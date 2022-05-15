#----------------------------------------------------------------------
# static components to plot summary results over all selected samples
#----------------------------------------------------------------------

# module ui function
analyzeSNVsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$analyzeSNVs)

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
        
        # sortable table of SVs
        fluidRow(
            filteredSvsTableUI(ns, width = 12)
        ),
        
        # expanded views of a single junction
        fluidRow(
            box(
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE,
                title = tagList(
                    "Junction Alignment to Source Alleles",
                    settingsUI(ns('alignmentSettings'))
                ),
                tags$style(HTML("
                    .junction { background-color: #ddd;}
                    .alignment .base_A { color: green; }
                    .alignment .base_C { color: blue; }
                    .alignment .base_G { color: brown; }
                    .alignment .base_T { color: red; }
                    .matches .mismatch { color: red; }
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
