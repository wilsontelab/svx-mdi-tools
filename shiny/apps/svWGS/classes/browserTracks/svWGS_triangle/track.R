#----------------------------------------------------------------------
# svWGS_triangle trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svWGS_triangleTrackBuffer <- reactiveValues()
svWGS_triangleExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svWGS_triangleTrack <- function(...) {
    new_svx_triangleTrack(..., expandReactive = svWGS_triangleExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svWGS_triangleTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.svWGS_triangleTrack <- function(...){
    build.svx_triangle_track(
        ..., 
        trackBuffer = svWGS_triangleTrackBuffer, 
        loadFn      = svWGS_loadJunctions
    )
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svWGS_triangleTrack <- function(...){
    svx_junctionNavTable(
        ..., 
        expandReactive  = svWGS_triangleExpand, 
        loadFn          = svWGS_loadJunctions, 
        navTableFn      = svWGS_navTable_display, 
        expandFn        = svWGS_expandJunction
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svWGS_triangleTrack <- function(...){
    svx_handleJunctionClick(
        ..., 
        buffer          = svWGS_triangleTrackBuffer, 
        expandReactive  = svWGS_triangleExpand, 
        expandFn        = svWGS_expandJunction, 
        summarizeFn     = svWGS_summarizeJunctions,
        distType        = "triangle"
    )   
}

# expand method for the S3 class
expand.svWGS_triangleTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, expandReactive = svWGS_triangleExpand)
}

# expand2 method for the S3 class
expand2.svWGS_triangleTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
