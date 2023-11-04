#----------------------------------------------------------------------
# svWGS_nodes trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svWGS_nodesTrackBuffer <- reactiveValues()
svWGS_nodesExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svWGS_nodesTrack <- function(...) {
    new_svx_nodesTrack(..., expandReactive = svWGS_nodesExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svWGS_nodesTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.svWGS_nodesTrack <- function(...){
    build.svx_nodes_track(
        ..., 
        trackBuffer = svWGS_nodesTrackBuffer, 
        loadFn      = svWGS_loadJunctions,
        idCol = "SV_ID"
    )
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svWGS_nodesTrack <- function(...){
    svx_junctionNavTable(
        ..., 
        expandReactive  = svWGS_nodesExpand, 
        loadFn          = svWGS_loadJunctions, 
        navTableFn      = svWGS_navTable_display, 
        expandFn        = svWGS_expandJunction
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svWGS_nodesTrack <- function(...){
    svx_handleJunctionClick(
        ..., 
        buffer          = svWGS_nodesTrackBuffer, 
        expandReactive  = svWGS_nodesExpand, 
        expandFn        = svWGS_expandJunction, 
        summarizeFn     = svWGS_summarizeJunctions,
        distType        = "nodes"
    )                         
}

# expand method for the S3 class
expand.svWGS_nodesTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, expandReactive = svWGS_nodesExpand)
}

# expand2 method for the S3 class
expand2.svWGS_nodesTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
