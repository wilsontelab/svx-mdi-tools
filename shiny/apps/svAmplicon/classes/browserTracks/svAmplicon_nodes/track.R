#----------------------------------------------------------------------
# svAmplicon_nodes trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svAmplicon_nodesTrackBuffer <- reactiveValues()
# svAmplicon_nodesExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svAmplicon_nodesTrack <- function(trackId) {
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
items.svAmplicon_nodesTrack <- function(...) svAmplicon_showAmpliconsDialog(...)

# build method for the S3 class; REQUIRED
build.svAmplicon_nodesTrack <- function(...){
    build.svx_nodes_track(..., svAmplicon_nodesTrackBuffer, svAmplicon_loadJunctions, idCol = "jxnUniqueKey", 
                          isMultiSample = FALSE, sampleNameFn = svAmplicon_legendSampleNames, jxnFilterFn = svAmplicon_filterJunctions)
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svAmplicon_nodesTrack <- function(...){
    svx_ampliconNavTable(...)
}

# # plot interaction methods for the S3 class
# # called by trackBrowser if track$click, $hover, or $brush is TRUE, above
# click.svAmplicon_nodesTrack <- function(...){
#     svx_handleJunctionClick(
#         ..., svAmplicon_nodesTrackBuffer, 
#         svAmplicon_nodesExpand, svAmplicon_expandJunction, svAmplicon_summarizeJunctions,
#         distType = "nodes"
#     )                           
# }

# # expand method for the S3 class
# expand.svAmplicon_nodesTrack <- function(track, reference, coord, layout){
#     svx_handleJunctionExpansion(track, layout, svAmplicon_nodesExpand)
# }

# # expand2 method for the S3 class
# expand2.svAmplicon_nodesTrack <- function(track, reference, coord, selectedRowData){
#     svx_handleJunctionExpansion2(track, reference, selectedRowData)
# }
