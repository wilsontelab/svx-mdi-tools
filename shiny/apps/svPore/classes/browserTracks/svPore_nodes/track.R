#----------------------------------------------------------------------
# svPore_nodes trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svPore_nodesTrackBuffer <- reactiveValues()
svPore_nodesExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svPore_nodesTrack <- function(...) {
    new_svx_nodesTrack(..., svPore_nodesExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_nodesTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.svPore_nodesTrack <- function(...){
    build.svx_nodes_track(..., svPore_nodesTrackBuffer, svPore_loadJunctions, idCol = "clusterN")
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svPore_nodesTrack <- function(...){
    svx_junctionNavTable(..., svPore_nodesExpand, svPore_loadJunctions, svPore_navTable_display, svPore_expandJunction)
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svPore_nodesTrack <- function(...){
    svx_handleJunctionClick(
        ..., svPore_nodesTrackBuffer, 
        svPore_nodesExpand, svPore_expandJunction, svPore_summarizeJunctions,
        distType = "nodes"
    )                           
}

# expand method for the S3 class
expand.svPore_nodesTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, svPore_nodesExpand)
}

# expand2 method for the S3 class
expand2.svPore_nodesTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
