#----------------------------------------------------------------------
# static components to show plots with QC values and distributions
# for all samples in a sampleSet, in order to explore data quality
#----------------------------------------------------------------------

# module ui function
libraryQCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)    

    # override missing options to defaults
    options <- setDefaultOptions(options, stepModuleInfo$libraryQC)

    # incorporate options text into templates
    leaderText <- tagList(
        tags$p(HTML(options$leaderText))
    )
    
    # unpad columns with boxes
    unpad <- "padding-left: 0; padding-right: 0;"

    # the available QC metrics
    libraryMetricTypes <- c(
        "maxTLen",
        "nReadPairs",
        "nSourceMolecules",
        "onTargetCoverage",
        "offTargetCoverage",
        "enrichment",
        "efficiency",
        "nSvs"
    )

    # return the UI contents
    standardSequentialTabItem(
        HTML(paste( options$longLabel, settingsUI(ns('settings')) )),
        # HTML(options$longLabel),
        leaderText,

        # box for selecting sample sources
        fluidRow( box(
            width = 12,
            column(
                width = 12,
                sampleSetUI(ns('sampleSet'))
            ),
            column(
                width = 12,
                # bsButton(ns('applyFailedFilters'), 'Apply Failure Settings',  style = "warning"), 
                # HTML('&nbsp;'),
                # bsButton(ns('clearFailedMarks'),   'Clear All Failure Marks', style = "success"),
                # HTML('&nbsp;'),
                textOutput(ns('nFailedLibraries'), inline = TRUE)           
            )
        ) ),
        
        # examination of quality metrics
        fluidRow(
            column(
                width = 6,
                style = unpad,
                box(width = 12,
                    title = "Library Metrics",
                    interactiveBarplotUI(ns('libraryMetricsPlot'), height = '400px')
                )
            ),
            column(
                width = 6,
                style = unpad,
                box(
                    width = 12,
                    title = "Metric Relationships",
                    collapsible = TRUE,
                    collapsed = FALSE,
                    fluidRow(
                        column(
                            width = 8,
                            interactiveScatterplotUI(ns('metricRelationshipsPlot'), height = '400px')
                        ),
                        column(
                            width = 4,
                            selectInput(ns('yAxisMetric'), 
                                        "Y Axis", 
                                        choices = libraryMetricTypes, 
                                        selected = "grouped"),     
                            tags$p("... as a function of ..."),                   
                            selectInput(ns('xAxisMetric'), 
                                        "X Axis", 
                                        choices = libraryMetricTypes, 
                                        selected = "alignRate")
                        )
                    )            
                )
            )
        ),
        
        # box with sortable table of samples with FAILED checkboxes
        fluidRow(
            column(
                width = 12,
                style = unpad,
                box(width = 12,
                    bufferedTableUI(ns('librariesTable'))
                )
            )
        ),

        # explore individual sample distributions
        fluidRow(
            column(
                width = 6,
                style = unpad,
                box(width = 12,
                    title = "Insert Sizes",
                    plotOutput(ns('insertSizesPlot'), height = '400px')
                )
            ),
            column(
                width = 6,
                style = unpad,
                box(width = 12,
                    title = "Strand Family Sizes",
                    plotOutput(ns('familySizesPlot'), height = '400px')
                )
            ) 
        )
    )    
}
