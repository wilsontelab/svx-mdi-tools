#----------------------------------------------------------------------
# UI components for the collapsibleDiv widget module
#----------------------------------------------------------------------

# module ui function
collapsibleDivUI <- function(id, ..., maxWidth = "700px", collapsed = TRUE) {
    ns <- NS(id)
    whiteSpace <- if(collapsed) "nowrap" else "normal"
    id <- ns("collapsibleDiv")
    cssId <- paste0("#", id) 
    tags$div(
        id = id,
        ...,
        style = paste("overflow: hidden; text-overflow: ellipsis; cursor: pointer; max-width:", maxWidth, "; white-space: ", whiteSpace, ";"),
        tags$script(paste0(
            "$('", 
            cssId, 
            "').on('click', function(){ $('", cssId, "').css('white-space', $('", cssId, "').css('white-space') == 'nowrap' ? 'normal' : 'nowrap') })"
        ))
    )
}
