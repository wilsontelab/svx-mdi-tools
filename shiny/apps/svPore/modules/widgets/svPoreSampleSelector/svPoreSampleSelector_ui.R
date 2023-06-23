#----------------------------------------------------------------------
# UI components for the svPoreSampleSelector widget module
#----------------------------------------------------------------------

# module ui function
svPoreSampleSelectorUI <- function(id) {
    ns <- NS(id)
    fluidRow(
        dataSourceTableUI(ns("source"), "Source", width = 8, collapsible = FALSE, inFluidRow = FALSE),
        bufferedTableUI(ns("samples"), title = "Sample(s)", width = 4, collapsible = FALSE,
                        status = "primary", solidHeader = TRUE)
    )
}
