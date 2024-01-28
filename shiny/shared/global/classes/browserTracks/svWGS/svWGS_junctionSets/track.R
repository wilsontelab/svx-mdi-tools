#----------------------------------------------------------------------
# svWGS_junctionSets trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
svWGS_junctionSetsTrackBuffer <- reactiveValues()
svWGS_junctionSetsExpand <- reactiveVal(NULL)

# constructor for the S3 class
new_svWGS_junctionSetsTrack <- function(...) {
    new_svx_navigatorTrack(..., expandReactive = svWGS_junctionSetsExpand)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.svWGS_junctionSetsTrack <- function(track, session, ...){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select Junction Set",
        itemTypePlural = "Junction Sets",
        tableData = reactive(app$junctionSets$getSavedJunctionSets()[, .(
            Name,
            Sample,
            Genome,
            SV_Types,
            N_Samples,
            N_Junctions,
            Hash
        )]),
        keyColumn = "Hash",
        extraColumns = c("Name", "Sample", "Genome", "SV_Types","N_Samples"),
        size = "xl"
    )
}

# build method for the S3 class; REQUIRED
build.svWGS_junctionSetsTrack <- function(...){
    NULL
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
navigation.svWGS_junctionSetsTrack <- function(track, ...){
    svx_junctionNavTable(
        track,
        ..., 
        expandReactive  = svWGS_junctionSetsExpand, 
        loadFn          = NULL, 
        navTableFn      = svWGS_navTable_display, 
        expandFn        = svWGS_expandJunction,
        jxnsReactive = reactive({
            junctionSets <- track$settings$items()
            req(junctionSets)
            app$junctionSets$getJunctionSetSvs(names(junctionSets))
        })
    )
}

# expand method for the S3 class
expand.svWGS_junctionSetsTrack <- function(track, reference, coord, layout){
    svx_handleJunctionExpansion(track, layout, expandReactive = svWGS_junctionSetsExpand)
}

# expand2 method for the S3 class
expand2.svWGS_junctionSetsTrack <- function(track, reference, coord, selectedRowData){
    svx_handleJunctionExpansion2(track, reference, selectedRowData)
}
