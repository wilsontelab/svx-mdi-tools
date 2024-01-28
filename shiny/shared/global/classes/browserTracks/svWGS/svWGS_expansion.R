# svWGS junction expansion

# output table definitions
svWGS_objectTable <- function(jxn){
    req(jxn)
    jxn[, .(
        svId = SV_ID,
        type = svx_jxnType_codeToX(edgeType, "name"), # edgeType converted to svx style code by svWGS_loadJunctions
        size = SV_SIZE,
        insertSize = -MICROHOM_LEN,
        mapQ = mapQ,
        nSamples,
        nInstances,
        nSequenced = N_SPLITS,
        flnkCN = round(outerFlankCN, 1),
        cnc    = round(flankCNC, 1),           
        jxnCN  = round(junctionCN, 1),
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
            molId = MOL_ID,
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
svWGS_expandJunction_ <- function(jxn, track = NULL, targetId = NULL,  # called by any appStep that can expand junctions
                                  objectTableFn = NULL, expansionTableFn = NULL, expansionUIFn = NULL,
                                  Pixels_Per_Base = 2, Bases_Per_Line = 100, Alignment_Mode = "Reference Molecule",
                                  session = NULL, expandedObjectFn = NULL){

    # collect molecules that provided evidence for the junction call
    if(is.null(targetId)) targetId <- jxn$targetId
    molecules <- svWGS_loadMolecules(targetId, jxn$SV_ID)

    # store the expanding junction for possible saveJunction actions
    if(!is.null(expandedObjectFn)) expandedObjectFn(list(
        targetId = targetId,
        jxn = jxn
    ))

    # write the one-line object table with additional junction metadata
    if(!is.null(objectTableFn)){
        jxn %>% 
        svWGS_getJunction() %>% 
        svWGS_objectTable() %>% 
        objectTableFn() # junction metadata   
    }

    # write the multi-line table with one line per supporting molecule
    if(!is.null(expansionTableFn)){
        molecules %>% 
        svWGS_expansionTable(jxn) %>% 
        expansionTableFn()
    }

    # if the junction was sequenced, create the expansion2 elements
    # a map and an alignment at base-level detail
    if(!is.null(expansionUIFn)){
        if(is.null(session)) session <- sessionSession
        junctionMap <- tryCatch({
            getJunctionMap(list(sv = jxn, mols = molecules[sample.int(.N)]))
        }, error = function(e) {
            stopSpinner(session)
            NULL
        })
        if(is.null(junctionMap)){
            expansionUIFn("")
        } else {
            startSpinner(session, message = "analyzing junction")
            expansionUIFn(tagList(
                if(svx_hasSavedJunctions()) fluidRow( 
                    style = "margin-bottom: 10px;",
                    column(
                        offset =4,
                        width = 4,
                        bsButton(
                            session$ns("saveJunction"), # handled by appServer >> observeEvent(input$saveJunction)
                            "Save Junction",
                            style = "primary",
                            block = TRUE
                        )
                    ),
                    column(
                        width = 4,
                        uiOutput(session$ns("saveJunctionFeedback"))
                    )
                ) else NULL,
                tryCatch({       
                    junctionMapTrackExpansionUI(junctionMap, Pixels_Per_Base) 
                }, error = function(e) ""),
                tryCatch({ 
                    junctionAlignmentTrackExpansionUI(junctionMap, Bases_Per_Line, Alignment_Mode) 
                }, error = function(e) { print(e); "" })
            ))     
        }
    }
}
svWGS_expandJunction <- function(jxn, track, layout){ # called from a browser track
    req(jxn)
    svWGS_expandJunction_(
        jxn = jxn,
        track = track,
        objectTableFn    = app$browser$objectTableData, 
        expansionTableFn = app$browser$expansionTableData, 
        expansionUIFn    = app$browser$expansionUI,
        Pixels_Per_Base  = getTrackSetting(track, "Junctions", "Pixels_Per_Base", 2),
        Bases_Per_Line   = getTrackSetting(track, "Junctions", "Bases_Per_Line", 100),
        Alignment_Mode   = getTrackSetting(track, "Junctions", "Alignment_Mode", "Reference Molecule"),
        expandedObjectFn = app$browser$expandedObjectData
    )
    stopSpinner(session)

    # # collect molecules that provided evidence for the junction call
    # molecules <- svWGS_loadMolecules(jxn$targetId, jxn$SV_ID)

    # # write the one-line object table with additional junction metadata
    # jxn %>% 
    # svWGS_getJunction() %>% 
    # svWGS_objectTable() %>% 
    # app$browser$objectTableData() # junction metadata    

    # # write the multi-line table with one line per supporting molecule
    # molecules %>% 
    # svWGS_expansionTable(jxn) %>% 
    # app$browser$expansionTableData()

    # # if the junction was sequence, create the expansion2 elements
    # # a map and an alignment at base-level detail
    # junctionMap <- tryCatch({
    #     getJunctionMap(list(sv = jxn, mols = molecules[sample.int(.N)]))
    # }, error = function(e) {
    #     stopSpinner(session)
    #     NULL
    # })
    # if(is.null(junctionMap)){
    #     app$browser$expansionUI("")
    # } else {
    #     startSpinner(session, message = "analyzing junction")
    #     app$browser$expansionUI(tagList(
    #         tryCatch({       junctionMapTrackExpansionUI(track, junctionMap) }, error = function(e) ""),
    #         tryCatch({ junctionAlignmentTrackExpansionUI(track, junctionMap) }, error = function(e) { print(e); "" })
    #     ))     
    # }

    req(FALSE) # do not show an in-browser track expansion, it is all too big and shown in expansionUI
}
