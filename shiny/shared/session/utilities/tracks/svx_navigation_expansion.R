# construct a junction cluster trackNav table
svx_junctionNavTable <- function(track, session, browserId, reference, coord, expandReactive){
    navTableName <- initTrackNav(track, session, "navTable") # table reactive functions are provided below  
    trackNavDataUnformatted <- reactive({
        selectedSources <- getSourcesFromTrackSamples(track$settings$items())
        jxns <- svx_getTrackJunctions(track, selectedSources)
        req(nrow(jxns) <= 10000)
        jxns
    })
    trackNavDataFormatted <- reactive({
        svx_navTable_display_app(trackNavDataUnformatted())
    })
    handleRowClick <- function(selectedRow){
        req(selectedRow)
        jxn <- trackNavDataUnformatted()[selectedRow]
        expandFn <- function(jobId){
            expandReactive(list(jxn = jxn, fnName = "svx_expandJunction_app"))
            app$browser$expandingTrackId(track$id)  
        }
        if(jxn$edgeType == svx_edgeTypes$TRANSLOCATION){
            handleTrackNavTableClick(track, "all", abs(jxn$node1), abs(jxn$node2), expandFn = expandFn)
        } else {
            handleTrackNavTableClick(track, jxn$cChrom1, jxn$cRefPos1, jxn$cRefPos2, expandFn = expandFn)
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

# appropriately disperse a nodes or triangle plot click event
svx_handleJunctionClick <- function(track, click, buffer, expandReactive, distFn = NULL, distType = NULL){
    jxns <- buffer[[track$id]]
    req(nrow(jxns) > 0) 
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
        expandReactive(list(jxn = jxn, fnName = "svx_expandJunction_app"))
        app$browser$expandingTrackId(track$id)  
    } else if(click$keys$shift){
        expandReactive(list(jxn = jxns, fnName = "svx_summarizeJunctions_app"))
        app$browser$expandingTrackId(track$id) 
    } else {
        svx_jumpToJunction(jxn)       
    }
}
svx_handleJunctionExpansion <- function(track, layout, expandReactive){
    x <- expandReactive()
    get(x$fnName)(x$jxn, track, layout)
}

# execute second-level expansion of a specific read, from an expansion table click
svx_handleJunctionExpansion2 <- function(track, reference, selectedRowData){
    app$browser$expansionUI(selectedRowData)
}
