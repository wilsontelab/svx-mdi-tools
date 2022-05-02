#----------------------------------------------------------------------
# static components to generate interactive plots for visualizing 
# batch effect normalization in bin count data 
#----------------------------------------------------------------------

# module ui function
normalizeBatchUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalizeBatch)

    # incorporate options text into templates
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )
    
    # unpad columns with boxes
    unpad <- "padding-left: 0; padding-right: 0;"
    failed <- "color: rgb(200,0,0); font-size: 1.2em; text-align: center;"
    
    # return the UI contents
    standardSequentialTabItem(
        #HTML(paste( options$longLabel, settingsUI(ns('settings')) )),
        options$longLabel,
        leaderText,
        
        # top box for selecting sample sources
        fluidRow( box(
            
            # select the Sample Set and Source (i.e., set of related samples/cells)
            width = 12,
            column(
                width = 4,
                sampleSetUI(ns('sampleSet'))
            ),
            column(
                width = 4,
                bookmarkInput('selectInput', ns('sourceId'), 'Sample Source', choices = c('-'))
            ),
            
            # select a specific cell from a set of related samples/cells
            column( 
                width = 4,
                bookmarkInput('listStepperButtonsUI', ns('stepper'))
            )
        ) ),
        
        fluidRow(
            
            # left column of boxes for QC summary plots
            column(width = 3,
                style = unpad,
                box(
                    width = 12,
                    title = "Copy Number Density",
                    interactiveScatterplotUI(ns('readDepthDensityPlot'), height = '250px')
                ),                 
                box(width = 12, 
                    title = "Sample vs. Median",
                    collapsible = TRUE,
                    collapsed = TRUE,
                    interactiveScatterplotUI(ns('sampleVsMedian'), height = '250px')
                )          
            ),
            
            # right column of boxed for genomic bin plots
            column(width = 9,
                style = unpad,
                box(
                    width = 12,
                    title = "Before (top) and after (bottom) batch effect normalization",
                    tagList(
                        div(style = failed, textOutput(ns('failedQC'))),
                        interactiveScatterplotUI(ns('gcCorrectedDepthPlot'), height = '250px'),
                        interactiveScatterplotUI(ns('beCorrectedDepthPlot'), height = '250px') 
                    )
                )                 
            )
        )
    )    
}
