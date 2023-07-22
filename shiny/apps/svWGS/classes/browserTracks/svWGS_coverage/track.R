#----------------------------------------------------------------------
# svWGS_coverage trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_svWGS_coverageTrack <- function(...) {
    new_svx_coverageTrack(...)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svWGS_coverageTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.svWGS_coverageTrack <- function(...){
    build.svx_coverageTrack(...)
}
