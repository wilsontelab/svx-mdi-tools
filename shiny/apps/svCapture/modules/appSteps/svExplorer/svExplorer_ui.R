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
        HTML(paste( options$longLabel, stepSettingsUI(ns('settings')) )),
        # HTML(options$longLabel),
        leaderText,

        # box for selecting sample sources
        fluidRow( box(
            width = 12,
            column(
                width = 12,
                sampleSelectorUI(ns('sampleSelector'))
            )
        ) ),
        
        # composite views of all filtered SVs
        fluidRow(
            column(
                width = 6,
                style = unpad,
                box(width = 12,
                    title = "SV Locations",
                    interactiveScatterplotUI(ns('svLocationsPlot'), height = '400px')
                )
            ),
            column(
                width = 6,
                style = unpad,
                box(width = 12,
                    title = "SV Properties",
                    interactiveScatterplotUI(ns('svPropertiesPlot'), height = '400px')
                )
            )
        ),
        
        # box with sortable table of SVs
        fluidRow(
            column(
                width = 12,
                style = unpad,
                box(width = 12,
                    bufferedTableUI(ns('svsTable'))
                )
            )
        )
        # ,

        # # explore individual sample distributions
        # fluidRow(
        #     column(
        #         width = 6,
        #         style = unpad,
        #         box(width = 12,
        #             title = "Insert Sizes",
        #             plotOutput(ns('insertSizesPlot'), height = '400px')
        #         )
        #     ),
        #     column(
        #         width = 6,
        #         style = unpad,
        #         box(width = 12,
        #             title = "Strand Family Sizes",
        #             plotOutput(ns('familySizesPlot'), height = '400px')
        #         )
        #     ) 
        # )
    )    
}
