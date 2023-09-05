# construct a junction cluster trackNav table
svx_junctionNavTable <- function(track, session, browserId, reference, coord, 
                                 expandReactive, loadFn, navTableFn, expandFn){
    navTableName <- initTrackNav(track, session, "navTable") # table reactive functions are provided below  
    trackNavDataUnformatted <- reactive({
        selectedSources <- getSourcesFromTrackSamples(track$settings$items())
        jxns <- svx_getTrackJunctions(track, selectedSources, loadFn, chromOnly = FALSE)
        req(nrow(jxns) <= 10000)
        jxns
    })
    trackNavDataFormatted <- reactive({
        navTableFn(trackNavDataUnformatted())
    })
    handleRowClick <- function(selectedRow){
        req(selectedRow)
        jxn <- trackNavDataUnformatted()[selectedRow]
        expandFn_ <- function(jobId){
            expandReactive(list(d = jxn, fn = expandFn))
            app$browser$expandingTrack(list(trackId = track$id, object = jxn) )
        }
        if(jxn$edgeType == svx_edgeTypes$TRANSLOCATION){
            handleTrackNavTableClick(track, "all", abs(jxn$node1), abs(jxn$node2), expandFn = expandFn_)
        } else {
            handleTrackNavTableClick(track, jxn$cChrom1, jxn$cRefPos1, jxn$cRefPos2, expandFn = expandFn_)
        }
    }
    tagList(
        trackNavTable(
            track, 
            session, 
            browserId,
            navTableName, # the name as provided by initTrackNav
            tableData = trackNavDataFormatted, # populate a table based on track settings, etc.
            actionFn = handleRowClick
        )
    )  
}

# handle SV junction-based browser navigation
svx_jumpToJunction <- function(jxn){
    if(jxn$edgeType == svx_edgeTypes$TRANSLOCATION){
        app$browser$jumpToCoordinates("all", abs(jxn$node1), abs(jxn$node2))
    } else {
        app$browser$jumpToCoordinates(jxn$cChrom1, jxn$cRefPos1, jxn$cRefPos2)
    }
}

# construct an amplicon trackNav table
svx_ampliconNavTable <- function(track, session, browserId, reference, coord, 
                                 parseFn = NULL, isMultiSample = FALSE){
    navTableName <- initTrackNav(track, session, "navTable") # table reactive functions are provided below  
    trackNavData <- reactive({ # simple nav table, no need for formatted/unformatted
        svx_getTrackAmpliconTargets(track, parseFn = parseFn, isMultiSample = isMultiSample)
    })
    handleRowClick <- function(selectedRow){
        req(selectedRow)
        amplicon <- trackNavData()[selectedRow]
        if(is.na(amplicon$end)){
            amplicon$end   <- amplicon$start + 500
            amplicon$start <- amplicon$start - 500
        }
        handleTrackNavTableClick(track, amplicon$chrom, amplicon$start, amplicon$end, expandFn = NULL)
    }
    tagList(
        trackNavTable(
            track, 
            session, 
            browserId,
            navTableName, # the name as provided by initTrackNav
            tableData = trackNavData, # populate a table based on track settings, etc.
            actionFn = handleRowClick
        )
    )  
}

# appropriately disperse a nodes or triangle plot click event
svx_handleJunctionClick <- function(track, click, buffer, 
                                    expandReactive, expandFn, summarizeFn,
                                    distFn = NULL, distType = NULL){
    jxns <- buffer[[track$id]]
    req(nrow(jxns) > 0) 
    startSpinner(session)
    dist <- if(!is.null(distFn)) distFn(jxns) else switch(
        distType,
        triangle = jxns[, sqrt((center - click$coord$x) ** 2 + (size - click$coord$y) ** 2)],
        nodes = {
            y2 <- jxns[, ((y - click$coord$y) / click$coord$y) ** 2]
            jxns[, pmin(
                sqrt(((pos1 - click$coord$x) / click$coord$x) ** 2 + y2),
                sqrt(((pos2 - click$coord$x) / click$coord$x) ** 2 + y2)
            )]
        }
    )
    jxn <- jxns[which.min(dist)]  
    if(click$keys$ctrl){
        expandReactive(list(d = jxn, fn = expandFn))
        app$browser$expandingTrack(list(trackId = track$id, object = jxn) )
    } else if(click$keys$shift){
        expandReactive(list(d = jxns, fn = summarizeFn))
        app$browser$expandingTrack(list(trackId = track$id, object = jxns) )
    } else {
        svx_jumpToJunction(jxn)       
    }
}
svx_handleJunctionExpansion <- function(track, layout, expandReactive){
    x <- expandReactive()
    x$fn(x$d, track, layout)
}

# execute second-level expansion of a specific read, from an expansion table click
svx_handleJunctionExpansion2 <- function(track, reference, selectedRowData){
    app$browser$expansionUI(selectedRowData)
}
