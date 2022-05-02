#----------------------------------------------------------------------
# static components for the DESeq2 analysis type results viewer
#----------------------------------------------------------------------

# inline module ui function
DESeq2UI <- function(id, options) {
    
    # initialize namespace
    ns <- NS(id)

    # always wrap in a div that acts as the removeUI css target
    tags$div(

        # the primary table of DESeq2 results by gene or feature
        fluidRow(
            resultsTableUI(
                ns('genesTable'),
                title = 'Results By Gene'  
            )
        ),
        
        # composite plots of all genes
        fluidRow(
            collapsibleBox(
                title = "Volcano Plot",
                interactiveScatterplotUI(ns('volcanoPlot'))
            ),
            collapsibleBox(
                title = "MA Plot",
                interactiveScatterplotUI(ns('MAPlot'))
            )
        ),
        
        # plot of single gene, all samples
        fluidRow(
            collapsibleBox(
                title = "Single Gene Plot",
                interactiveScatterplotUI(ns('genePlot'))
            )
        )  
    )
}
