#----------------------------------------------------------------------
# UI components for svWGS explorer data selectors
#----------------------------------------------------------------------

# module ui function
svWGS_explorer_dataSelectorsUI <- function(id) {
    ns <- NS(id)
    fluidRow(
        dataSourceTableUI(
            ns("sources"), 
            "Data Source", 
            width = 6, 
            collapsible = TRUE,
            inFluidRow = FALSE
        ),
        bufferedTableUI(
            ns("sampleGenomes"),
            "Sample with Genome",
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE
        )
    )   
}
