# svWGS junction expansion

# output table definitions
svWGS_objectTable <- function(jxn){
    req(jxn)
    jxn[, .(
        svId = SV_ID,
        type = svx_jxnTypes[edgeType, name],
        size = SV_SIZE,
        insertSize = -MICROHOM_LEN,
        # mapQ,
        # identity = gapCompressedIdentity,
        # alnBaseQual = as.integer(alnBaseQual),
        nSamples,
        nInstances,
        rStart = paste0(cChrom1, ":", cRefPos1, if(cStrand1 == 1) "+" else "-"),
        rEnd   = paste0(cChrom2, ":", cRefPos2, if(cStrand2 == 1) "+" else "-")
    )]
}
svWGS_expansionTable <- function(edges){
    req(edges)
    edges[, .(
        readKey,
        channel,
        nStrands,
        edge = paste(edgeN, edgeType, sep = ":"),
        eventSize,
        insertSize,
        qStart,
        qEnd,
        rStart = paste0(chrom1, ":", refPos1, strand1),
        rEnd   = paste0(chrom2, ":", refPos2, strand2)
    )]
}

# show a detailed plot and table of the molecule support for a junction cluster
svx_expandJunction_app <- function(jxn, track, layout){
    req(jxn)

    padding <- padding(track, layout)
    height <- getBrowserTrackSetting(track, "Junctions", "Plot_Height", 3) # height in inches
    locusPadding <- getBrowserTrackSetting(track, "Junctions", "Locus_Padding", 50000)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svWGS_triangle zoom", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))

        jxn %>% 
        svWGS_getJunction() %>% 
        svWGS_objectTable() %>% 
        app$browser$objectTableData() # junction metadata

        plot(1:10)

        # edges <- loadClusterEdges(jxn$sourceId, jxn$clusterN) 

        # jxn %>% 
        # svWGS_getJunction() %>% 
        # svWGS_objectTable() %>% 
        # app$browser$objectTableData() # junction cluster metadata

        # edges %>% 
        # svWGS_expansionTable() %>% 
        # app$browser$expansionTableData() # supporting junctions metadata
        
        # edges %>% 
        # setAlignmentLoci(locusPadding) %>% 
        # parseEdgeForMolPlot() %>% 
        # renderMoleculePlot(height, getBrowserTrackSetting(track, "Junctions", "Plot_Speed", "fast")) # supporting junction image
    })
    stopSpinner(session)  

    # return the track's magick image and associated metadata
    list(
        # xlim  = xlim, # will be rescaled to c(0, 1) across all stacked zoom images
        ylim  = NA,
        mai   = mai,
        image = image
    )
}
