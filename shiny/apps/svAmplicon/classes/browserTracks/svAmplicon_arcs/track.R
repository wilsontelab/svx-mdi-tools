#----------------------------------------------------------------------
# svAmplicon_arcs trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svAmplicon_arcsTrackBuffer <- reactiveValues()
# svAmplicon_arcsExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svAmplicon_arcsTrack <- function(trackId) {
    list(
        click  = FALSE,
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
items.svAmplicon_arcsTrack <- function(...) svAmplicon_showAmpliconsDialog(...)

# build method for the S3 class; REQUIRED
build.svAmplicon_arcsTrack <- function(...){
    build.svx_arcs_track(..., svAmplicon_arcsTrackBuffer, svAmplicon_loadJunctions, idCol = "jxnUniqueKey", 
                         isMultiSample = FALSE, sampleNameFn = svAmplicon_legendSampleNames)
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svAmplicon_arcsTrack <- function(...){
    svx_ampliconNavTable(...)
}

# # plot interaction methods for the S3 class
# # called by trackBrowser if track$click, $hover, or $brush is TRUE, above
# click.svAmplicon_arcsTrack <- function(...){
#     svx_handleJunctionClick(
#         ..., svAmplicon_arcsTrackBuffer, 
#         svAmplicon_arcsExpand, svAmplicon_expandJunction, svAmplicon_summarizeJunctions,
#         distType = "nodes"
#     )                           
# }

# # expand method for the S3 class
# expand.svAmplicon_arcsTrack <- function(track, reference, coord, layout){
#     svx_handleJunctionExpansion(track, layout, svAmplicon_arcsExpand)
# }

# # expand2 method for the S3 class
# expand2.svAmplicon_arcsTrack <- function(track, reference, coord, selectedRowData){
#     svx_handleJunctionExpansion2(track, reference, selectedRowData)
# }
