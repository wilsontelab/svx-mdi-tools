# support functions for svPore track assembly

# selecting the samples to plot
svPore_trackItems <- function(track, session, input, reference){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select Samples",
        itemTypePlural = "Samples",
        tableData = reactive({
            uploadName <- appStepNamesByType$upload
            x <- as.data.table(app[[uploadName]]$outcomes$samples())
            req(x)
            x[, .(Sample_ID, Project)]
        }),
        keyColumn = "Sample_ID",
        extraColumns = c("Project"),
        # ,
        # options = list(
        #     XXX = list(
        #         type = "selectInput", # or textInput, etc.
        #         args = list(
        #             choices = c("aaa", "bbb"),
        #             selected = "aaa",
        #             width = "50px"                  
        #         )
        #     )
        # ),
        size = "l" # xl
    )
}

# show several plots that summarize the junction clusters of one or more samples
svPore_summarizeJunctionClusters <- function(jc, track, layout){
    req(jc)
    jc[, alnBaseQual := as.integer(alnBaseQual)]
    jc[, alnSize_kb  := as.integer(alnSize / 1000)]

    plotH <- function(dt, lwd = 2){
        plot(
            dt, 
            typ = "h", 
            ylim = c(0, max(dt[[2]]) * 1.05), 
            xlab = names(dt)[1],
            ylab = "N",
            lwd = lwd
        )
    }    

    padding <- padding(track, layout)
    height <- getBrowserTrackSetting(track, "Summarize", "Plot_Height", 3) # height in inches
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svPore_triangle summary", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        layout(matrix(1:6, nrow = 2, byrow = TRUE)) 
        par(mar = c(4.6, 4.1, 0.6, 0.1), cex = 1)

        # mapQ            = max(mapQ), # these values aggregate over all junction in the cluster, even the fuzzy matched ones
        # gapCompressedIdentity = max(gapCompressedIdentity),

        plotH(jc[, .N, keyby = .(nSamples)])
        plotH(jc[, .N, keyby = .(nInstances)]) # nInstances <= 50
        plotH(jc[between(insertSize, -50, 50), .N, keyby = .(insertSize)])
        plotH(jc[, .N, keyby = .(alnBaseQual)])
        plotH(jc[, .N, keyby = .(alnSize_kb)])
        barplot(N ~ edgeType, jc[, .N, keyby = .(edgeType)])
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

# show a detailed plot and table of the molecule support for a junction cluster
svPore_objectTable <- function(jc){
    req(jc)
    jc[, .(
        clusterN,
        edgeType,
        eventSize,
        insertSize,
        mapQ,
        identity = gapCompressedIdentity,
        baseQual,
        nSamples,
        nInstances,
        rStart = paste0(cChrom1, ":", cRefPos1, if(cStrand1 == 1) "+" else "-"),
        rEnd   = paste0(cChrom2, ":", cRefPos2, if(cStrand2 == 1) "+" else "-")
    )]
}
svPore_expansionTable <- function(edges){
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
svPore_expandJunctionCluster <- function(jc, track, layout){
    req(jc)
    padding <- padding(track, layout)
    height <- getBrowserTrackSetting(track, "Junctions", "Plot_Height", 3) # height in inches
    locusPadding <- getBrowserTrackSetting(track, "Junctions", "Locus_Padding", 50000)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svPore_triangle zoom", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))

        edges <- loadClusterEdges(jc$sourceId, jc$clusterN) 

        jc %>% 
        getJunctionCluster() %>% 
        svPore_objectTable() %>% 
        app$browser$objectTableData() # junction cluster metadata

        edges %>% 
        svPore_expansionTable() %>% 
        app$browser$expansionTableData() # supporting junctions metadata
        
        edges %>% 
        setAlignmentLoci(locusPadding) %>% 
        parseEdgeForMolPlot() %>% 
        renderMoleculePlot(height, getBrowserTrackSetting(track, "Junctions", "Plot_Speed", "fast")) # supporting junction image
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

# handle SV junction based browser navigation
jumpToJunctionCluster <- function(jc){
    if(jc$edgeType == edgeTypes$TRANSLOCATION){
        app$browser$jumpToCoordinates("all", abs(jc$node1), abs(jc$node2))
    } else {
        app$browser$jumpToCoordinates(jc$cChrom1, jc$cRefPos1, jc$cRefPos2)
    }
}

# appropriately disperse a nodes or triangle click event
handleJunctionClustersClick <- function(track, click, buffer, expandReactive, distFn){
    jcs <- buffer[[track$id]]
    req(nrow(jcs) > 0)  
    dist <- distFn(jcs)
    jc <- jcs[which.min(dist)]  
    if(click$keys$ctrl){
        expandReactive(list(jc = jc, fnName = "svPore_expandJunctionCluster"))
        app$browser$expandingTrackId(track$id)  
    } else if(click$keys$shift){
        expandReactive(list(jc = jcs, fnName = "svPore_summarizeJunctionClusters"))
        app$browser$expandingTrackId(track$id) 
    } else {
        jumpToJunctionCluster(jc)       
    }
}
handleJunctionClusterExpansion <- function(track, layout, expandReactive){
    x <- expandReactive()
    get(x$fnName)(x$jc, track, layout)
}
