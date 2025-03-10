#----------------------------------------------------------------------
# amplicon_nodes trackBrowser track (i.e., a browserTrack)
# same structure and function as svWGS_junctions but now filtered to only show junctions
# within or between amplicon_highCopyRegions
#----------------------------------------------------------------------
amplicon_nodesTrackBuffer <- reactiveValues()
amplicon_nodesExpand <- reactiveVal(NULL)
# amplicon_loadJunctions(targetId){
#     jxns <- svWGS_loadJunctions(targetId)
#     # likely use data.table intersect function here
#     # denote junctions as within or between (also, 1-end in?)
#     jxns
# }

# constructor for the S3 class
new_amplicon_nodesTrack <- function(...) {
    new_svx_nodesTrack(..., expandReactive = amplicon_nodesExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.amplicon_nodesTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.amplicon_nodesTrack <- function(...){
    build.svx_nodes_track(
        ..., 
        trackBuffer = amplicon_nodesTrackBuffer, 
        loadFn      = amplicon_loadJunctions,
        idCol = "SV_ID"
    )
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.amplicon_nodesTrack <- function(...){
    svx_junctionNavTable(
        ..., 
        expandReactive  = amplicon_nodesExpand, 
        loadFn          = amplicon_loadJunctions, 
        navTableFn      = svWGS_navTable_display, 
        expandFn        = svWGS_expandJunction
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.amplicon_nodesTrack <- function(...){
    svx_handleJunctionClick(
        ..., 
        buffer          = amplicon_nodesTrackBuffer, 
        expandReactive  = amplicon_nodesExpand, 
        expandFn        = svWGS_expandJunction, 
        summarizeFn     = svWGS_summarizeJunctions,
        distType        = "nodes"
    )                         
}

# expand method for the S3 class
expand.amplicon_nodesTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, expandReactive = amplicon_nodesExpand)
}

# expand2 method for the S3 class
expand2.amplicon_nodesTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
