# svWGS junction expansion

# output table definitions
svWGS_objectTable <- function(jxn){
    req(jxn)
    jxn[, .(
        svId = SV_ID,
        type = svx_jxnType_codeToX(edgeType, "name"),
        size = SV_SIZE,
        insertSize = -MICROHOM_LEN,
        mapQ = mapQ,
        nSamples,
        nInstances,
        nSequenced = N_SPLITS,
        flankLen = min(FLANK_LEN1, FLANK_LEN2),
        nLinkedJxns = nLinkedJunctions,
        rStart = paste0(cChrom1, ":", cRefPos1, if(cStrand1 == 1) "+" else "-"),
        rEnd   = paste0(cChrom2, ":", cRefPos2, if(cStrand2 == 1) "+" else "-")
    )]
}
svWGS_expansionTable <- function(molecules, jxn){
    req(molecules)
    cbind(
        molecules[, .(
            sample = SAMPLE,
            isRef = as.logical(IS_REFERENCE),            
            merge = IS_MERGED,
            nodeClass = tolower(svx_nodeClasses_rev[NODE_CLASS + 1]),
            strand = ifelse(MOL_STRAND == 1, "-", "+"),            
            mapQ1 = MAPQ_1,
            mapQ2 = MAPQ_2,
            innerNode1 = NODE_1,
            innerNode2 = NODE_2,
            outerPos1 = OUT_POS1,
            outerPos2 = OUT_POS2
        )]
    )
}

# show a detailed plot and table of the molecule support for a junction cluster
svWGS_expandJunction <- function(jxn, track, layout){
    req(jxn)

    # collect molecules that provided evidence for the junction call
    molecules <- svWGS_loadMolecules(jxn$sourceId, jxn$SV_ID)

    # write the one-line object table with additional junction metadata
    jxn %>% 
    svWGS_getJunction() %>% 
    svWGS_objectTable() %>% 
    app$browser$objectTableData() # junction metadata    

    # write the multi-line table with one line per supporting molecule
    molecules %>% 
    svWGS_expansionTable(jxn) %>% 
    app$browser$expansionTableData()

    # if the junction was sequence, create the expansion2 elements
    # a map and an alignment at base-level detail
    junctionMap <- tryCatch({
        getJunctionMap(list(sv = jxn, mols = molecules[sample.int(.N)]))
    }, error = function(e) {
        stopSpinner(session)
        NULL
    })
    if(is.null(junctionMap)){
        app$browser$expansionUI("")
    } else {
        startSpinner(session, message = "analyzing junction")
        app$browser$expansionUI(tagList(
            tryCatch({       junctionMapTrackExpansionUI(track, junctionMap) }, error = function(e) ""),
            tryCatch({ junctionAlignmentTrackExpansionUI(track, junctionMap) }, error = function(e) { print(e); "" })
        ))     
    }

    # do not show an in-browser track expansion, it is all too big and shown in expansionUI
    # leave code block below in case we come up with some else to plot as a track...
    req(FALSE)

    # # set the expansion track layout
    # padding <- padding(track, layout)
    # height <- getBrowserTrackSetting(track, "Junctions", "Junction_Plot_Height", 3) # height in inches

    # # use the mdiTrackImage helper function to create the track image
    # mai <- NULL
    # image <- mdiTrackImage(layout, height, message = "svWGS_triangle zoom", function(...){
    #     mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
    # })
    # stopSpinner(session)  

    # # return the track's magick image and associated metadata
    # list(ylim  = NA, mai   = mai, image = image)
}
