#----------------------------------------------------------------------
# static components to filter and examine SV locations and junction sequences
#----------------------------------------------------------------------

# module ui function
svWGSExplorerUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$svWGSExplorer)

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
                width = 8,
                title = "SV Locations"
            ),
            svPropertiesPlotUI(ns, width = 4)
        ),
        
        # sortable table of SVs
        fluidRow(
            filteredSvsTableUI(ns, width = 12)
        ),
        
        # expanded views of a single junction
        fluidRow(
            junctionMapUI(ns, width = 12)
        ),
        fluidRow(
            junctionNodesPlotUI(ns, width = 4),
            junctionAlignmentUI(ns, width = 8)          
        ),
        ""
    )    
}
