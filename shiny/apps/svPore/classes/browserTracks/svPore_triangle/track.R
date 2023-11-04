#----------------------------------------------------------------------
# svPore_triangle trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svPore_triangleTrackBuffer <- reactiveValues()
svPore_triangleExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svPore_triangleTrack <- function(...) {
    new_svx_triangleTrack(..., expandReactive = svPore_triangleExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svPore_triangleTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.svPore_triangleTrack <- function(...){
    build.svx_triangle_track(
        ..., 
        trackBuffer = svPore_triangleTrackBuffer, 
        loadFn      = svPore_loadJunctions
    )
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svPore_triangleTrack <- function(...){
    svx_junctionNavTable(
        ..., 
        expandReactive  = svPore_triangleExpand, 
        loadFn          = svPore_loadJunctions, 
        navTableFn      = svPore_navTable_display, 
        expandFn        = svPore_expandJunction
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
click.svPore_triangleTrack <- function(...){
    svx_handleJunctionClick(
        ..., 
        buffer          = svPore_triangleTrackBuffer, 
        expandReactive  = svPore_triangleExpand, 
        expandFn        = svPore_expandJunction, 
        summarizeFn     = svPore_summarizeJunctions,
        distType        = "triangle"
    )    
}

# expand method for the S3 class
expand.svPore_triangleTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, expandReactive = svPore_triangleExpand)
}

# expand2 method for the S3 class
expand2.svPore_triangleTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
