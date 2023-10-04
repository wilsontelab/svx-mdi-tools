#----------------------------------------------------------------------
# svAmplicon_amplicons trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svAmplicon_ampliconsTrackBuffer <- reactiveValues()
# svAmplicon_ampliconsExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svAmplicon_ampliconsTrack <- function(trackId) {
    list(
        click  = TRUE,
        hover  = FALSE,
        brush  = FALSE,
        items  = TRUE,
        expand = NULL, # expandReactive,
        expand2 = FALSE,
        navigation = TRUE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svAmplicon_ampliconsTrack <- function(...) svAmplicon_showAmpliconsDialog(...)

# build method for the S3 class; REQUIRED
build.svAmplicon_ampliconsTrack <- function(...){
    build.svx_amplicons_track(..., svAmplicon_ampliconsTrackBuffer, isMultiSample = FALSE)
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svAmplicon_ampliconsTrack <- function(...){
    svx_ampliconNavTable(..., isMultiSample = FALSE)
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svAmplicon_ampliconsTrack <- function(...){

    # svx_handleJunctionClick(
    #     ..., svAmplicon_ampliconsTrackBuffer, 
    #     svAmplicon_ampliconsExpand, svAmplicon_expandJunction, svAmplicon_summarizeJunctions,
    #     distType = "nodes"
    # )                           
}
