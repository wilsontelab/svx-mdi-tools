# svPore junction expansion

# output table definitions
svPore_objectTable <- function(jxn){
    req(jxn)
    jxn[, .(
        cluster = clusterN,
        type = svx_jxnType_codeToX(edgeType, "name"),
        size = eventSize,
        insertSize,
        mapQ,
        identity = gapCompressedIdentity,
        alnBaseQual = as.integer(alnBaseQual),
        nSamples,
        nInstances,
        rStart = paste0(cChrom1, ":", cRefPos1, if(cStrand1 == 1) "+" else "-"),
        rEnd   = paste0(cChrom2, ":", cRefPos2, if(cStrand2 == 1) "+" else "-")
    )]
}
svPore_expansionTable <- function(edges){
    req(edges)
    edges[, .(
        sample,
        channel,        
        readI,
        duplex,
        duplexI = duplexCluster,
        duplex2, 
        duplex2I = duplexCluster2,
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
svPore_expandJunction <- function(jxn, track, layout, ...){
    req(jxn)
    padding <- padding(track, layout)
    height <- getBrowserTrackSetting(track, "Junctions", "Junction_Plot_Height", 3) # height in inches
    locusPadding <- getBrowserTrackSetting(track, "Junctions", "Locus_Padding", 50000)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svPore_triangle zoom", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))

        edges <- svPore_loadEdges(jxn$targetId, jxn$clusterN) 

        jxn %>% 
        svPore_getJunction() %>% 
        svPore_objectTable() %>% 
        app$browser$objectTableData() # junction metadata

        edges %>% 
        svPore_expansionTable() %>% 
        app$browser$expansionTableData() # supporting edges (i.e., molecules) metadata
        
        edges %>% 
        setAlignmentLoci(locusPadding) %>% 
        parseEdgeForMolPlot() %>% 
        renderMoleculePlot(height, getBrowserTrackSetting(track, "Junctions", "Plot_Speed", "fast")) # supporting edges (i.e., molecules) image
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
