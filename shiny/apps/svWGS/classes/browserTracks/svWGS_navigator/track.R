#----------------------------------------------------------------------
# svWGS_navigator trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svWGS_navigatorTrackBuffer <- reactiveValues()
svWGS_navigatorExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svWGS_navigatorTrack <- function(...) {
    new_svx_navigatorTrack(..., expandReactive = svWGS_navigatorExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svWGS_navigatorTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.svWGS_navigatorTrack <- function(...){
    NULL
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svWGS_navigatorTrack <- function(...){
    svx_junctionNavTable(
        ..., 
        expandReactive  = svWGS_navigatorExpand, 
        loadFn          = svWGS_loadJunctions, 
        navTableFn      = svWGS_navTable_display, 
        expandFn        = svWGS_expandJunction
    )
}

# expand method for the S3 class
expand.svWGS_navigatorTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, expandReactive = svWGS_navigatorExpand)
}

# expand2 method for the S3 class
expand2.svWGS_navigatorTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
