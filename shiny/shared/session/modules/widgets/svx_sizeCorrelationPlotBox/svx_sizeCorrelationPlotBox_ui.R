#----------------------------------------------------------------------
# UI components for plotting SV insert size vs. SV size
#----------------------------------------------------------------------

# module ui function
svx_sizeCorrelationPlotBoxUI <- function(id, ...) {
    # ns <- NS(id)
    staticPlotBoxUI(
        id, 
        title = "SV Size vs. Insert Size",
        ...       
    )
}
