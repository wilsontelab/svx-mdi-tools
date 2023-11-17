#----------------------------------------------------------------------
# svx_navigator is a track type for creating a junctions navTable 
# without making any knind of track image
# in this way, you can track a subsets of SVs while visualizing all of them, etc.
#----------------------------------------------------------------------

# constructor
new_svx_navigatorTrack <- function(trackId, expandReactive) {
    list(
        click  = FALSE,
        hover  = FALSE,
        brush  = FALSE,
        items  = TRUE,
        expand = expandReactive,
        expand2 = FALSE,
        navigation = TRUE
    )
}

# build method for the S3 class; REQUIRED
build.svx_navigator_track <- function(...){
    NULL
}
