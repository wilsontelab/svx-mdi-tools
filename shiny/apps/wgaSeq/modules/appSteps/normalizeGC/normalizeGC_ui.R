#----------------------------------------------------------------------
# static components to generate interactive plots for optimizing GC bias
# normalization in bin count data and to reject low quality samples/cells
#----------------------------------------------------------------------

# module ui function
normalizeGCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalizeGC)

    # incorporate options text into templates
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )
    
    # unpad columns with boxes
    unpad  <- "padding-left: 0; padding-right: 0;"
    failed <- "color: rgb(200,0,0); font-size: 1.2em; text-align: center;"
    
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
            
            # inputs for user to assign attributes to each sample/cell
            column(
                width = 1,
                style = unpad, # best way to make an inline label for the input
                tags$strong(
                    style = "float:right; line-height: 2em;",
                    "Modal CN"
                )
            ),
            column( # the modal CN relevant to the set of selected bins, used to adjust output CN
                width = 1,
                style = "padding-left: 10px; padding-right: 0;",
                numericInput(ns('sampleModalCN'), label = NULL, value = 2, min = 1, max = 4, step = 1)
            ),        
            column( # flag that a sample is of unacceptably low quality (and would confound batch effect normalization)
                width = 3,
                div(
                    style = "display: inline-block; text-align: center; width: 100%; padding-top: 5px;",
                    radioButtons(ns('batchRejectSample'), label = NULL, inline = TRUE,
                                 choices = c("Batch"    = CONSTANTS$batched,
                                           "No_Batch"   = CONSTANTS$unbatched,
                                           "Reject   "  = CONSTANTS$rejected))
                )
            ),
            
            # select a specific cell from a set of related samples/cells
            column( 
                width = 7,
                bookmarkInput('listStepperButtonsUI', ns('stepper'))
            )
        ) ),
        
        # row of plots
        fluidRow(
            
            # left column of boxes for QC summary plots
            column(width = 3,
                style = unpad,
                box(width = 12,
                    collapsible = TRUE,
                    collapsed = FALSE,
                    title = "Copy Number Density",
                    tagList(
                        div(style = failed, textOutput(ns('failedQC1'))),
                        interactiveScatterplotUI(ns('gcReadDepthDensityPlot'), height = '250px')
                    )
                ),                     
                box(width = 12,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    title = "Bin Read Depth vs. GC",
                    interactiveScatterplotUI(ns('readDepthByGCPlot'), height = '250px')
                ),
                box(width = 12,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    title = "Sample Read Depths",
                    interactiveScatterplotUI(ns('sampleReadDepthsPlot'), height = '250px')
                ),   
                box(width = 12,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    title = "Lorenz Plot, All Samples",
                    interactiveScatterplotUI(ns('lorenzPlot'), height = '250px')
                )                       
            ),
            
            # right column of boxes for genomic bin plots
            column(width = 9,
                style = unpad,
                box(width = 12,
                    collapsible = TRUE,
                    collapsed = FALSE,
                    title = HTML(paste(tags$strong("Swipe in top plot"),
                                       "to select regions that have the Modal CN;",
                                       tags$strong("click in bottom plot"),
                                       "to mark CNVs")),
                    tagList(
                        div(style = failed, textOutput(ns('failedQC2'))),
                        interactiveScatterplotUI(ns('gc0ReadDepthPlot'), height = '250px')
                    ),          
                    interactiveScatterplotUI(ns('gcReadDepthPlot'),  height = '250px') 
                )                    
            )
        ),
        
        # row for the table of all samples in a SampleSet
        fluidRow(
            column(
                width = 12,
                style = unpad,
                box(width = 12,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    solidHeader = TRUE,
                    title = "Summary of all samples in sample set/group/type",
                    bufferedTableUI(ns('samplesTable'))
                )
            )
        )
    )    
}
