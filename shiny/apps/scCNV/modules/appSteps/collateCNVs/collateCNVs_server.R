#----------------------------------------------------------------------
# server components for the collateCNVs appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
collateCNVsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'collateCNVs'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
cnvGroupAssignments <- reactiveVal({ # this module's outcomes, assignment of CNVs into groups presumed identical (by descent?)
    dt <- data.table(
        groupKey = character(), 
        pendingKey = character()
    )
    dt[, key := character()]
    list(
        table = dt,
        entropy = sample(1e8, 1)
    )
})

#----------------------------------------------------------------------
# load source data
#----------------------------------------------------------------------
sourceIds <- dataSourceTableServer(
    "source", 
    selection = "multiple"
)
observeEvent(sourceIds(), {
    updateSampleSelector(session, input, sourceIds())
})
nSelectedSources <- reactive({ getNSelectedSources(input, sourceIds) })
cnvsTableData <- reactive({ # a table of CNVs that were kept from the selected sources
    getFilteredCnvData(input, sourceIds, settings)
})
cnvTableWithGroups <- reactive({ # combine cnvs and group assignments into one working table
    cnvs <- cnvsTableData()
    req(cnvs)
    merge(cnvs, cnvGroupAssignments()$table, by = "key", all.x = TRUE)
})

#----------------------------------------------------------------------
# create a display table of all kept CNVs across all selected data sources
#----------------------------------------------------------------------
cnvsTableDataPretty <- reactive({
    cnvs <- cnvsTableData()[, .SD, .SDcols = c(
        "project",
        "sampleName",
        "cell_id",
        "chrom",
        "start",
        "end",
        "size",
        'aneuploid',
        'telomeric',
        "cellCN",
        "referenceCN",
        "nWindows",
        "windowPower"
    )]
})
cnvsTable <- bufferedTableServer(
    "cnvsTable",
    id,
    input,
    tableData = cnvsTableDataPretty,
    selection = 'single',
    options = list()
)
observeEvent(cnvsTable$rows_selected(), {
    i <- cnvsTable$rows_selected()
    req(i)
    cnvs <- cnvsTableData()
    zoomChrom( cnvs[i, chrom] )
})

#----------------------------------------------------------------------
# plot all cells with CNVs, on one chromosome, with clickable areas for grouping
#----------------------------------------------------------------------
# invalidateZoomPlots <- reactiveVal(0)
invalidateZoomedCells <- character()    
zoomChrom <- reactiveVal("chr1")
output$chromPlots <- renderUI({ 
    cnvs <- cnvTableWithGroups()
    req(cnvs)
    zoomChrom <- zoomChrom()   
    req(zoomChrom)
    startSpinner(session, message = "loading source CNVs")
    sourceIds <- unique(cnvs$sourceId)
    projects <- lapply(sourceIds, getScCnvProjectData)
    names(projects) <- sourceIds
    startSpinner(session, message = "loading chrom plots")
    labelRow <- createZoomLabelRow(session, zoomChrom, commit = TRUE, leftPad = TRUE)
    cnvs <- cnvs[cnvs$chrom == zoomChrom]
    nCnvs <- nrow(cnvs)
    cellIds <- sapply(sourceIds, function(sourceId) list())
    nCells <- 0
    filledCells <- integer()
    cells <- if(nCnvs == 0) list() else lapply(1:nCnvs, function(i){
        cnv <- cnvs[i]  
        if(is.null(cellIds[[cnv$sourceId]][[cnv$cell_id]])){                   
            project <- projects[[cnv$sourceId]]
            cell <- project$cells[[cnv$cell_id]]
            cellIds[[cnv$sourceId]][[cnv$cell_id]] <<- TRUE 
            nCells <<- nCells + 1   
            filledCells <<- c(filledCells, i)         
            plotOneCellUI_match(cnv$sourceId, project, cell, settings, zoomChrom, invalidateZoomedCells, 
                                isMatchedCnv = isMatchedCnv, cnvs[sourceId == cnv$sourceId & cell_id == cnv$cell_id])
        } else NULL
    })
    cells <- cells[filledCells]
    order <- if(nCells > 2){
        I <- sapply(sourceIds, function(sourceId) length(cellIds[[sourceId]]) > 0)
        wps <- sapply(sourceIds[I], function(sourceId) projects[[sourceId]]$colData[bad == FALSE, min(windowPower)])
        minWp <- min(wps) 
        matrices <- lapply(sourceIds[I], getScCnvCellMatrix, settings, minWp)
        names(matrices) <- sourceIds[I]
        chromWindows <- projects[[1]]$windows[[minWp + 1]][, chrom == zoomChrom]
        cnMatrix <- do.call(cbind, lapply(sourceIds[I], function(sourceId) {
            matrices[[sourceId]][chromWindows, names(cellIds[[sourceId]])] 
        }))
        hclust(dist(t(cnMatrix)))$order
    } else 1:nCells
    isolate({ 
        nWorkingCells( nCells )
        initGenomePlotClicks(initGenomePlotClicks() + 1) 
    })
    invalidateZoomedCells <<- character()    
    stopSpinner(session)
    tagList(labelRow, cells[order])
})
nWorkingCells <- reactiveVal(0)

#----------------------------------------------------------------------
# handle actions to navigate through chromosomes
#----------------------------------------------------------------------
prevNextZoomChrom <- function(inc){
    cnvs <- cnvsTableData()
    zoomChrom <- zoomChrom()    
    chroms <- cnvs[, unique(chrom)]
    i <- which(chroms == zoomChrom) + inc
    if(length(i) == 0) i <- 1
    if(i < 1 || i > length(chroms)) return(NULL)
    zoomChrom(chroms[i])
}
observeEvent(input$prevZoomChrom, { prevNextZoomChrom(-1) })
observeEvent(input$nextZoomChrom, { prevNextZoomChrom( 1) })
observeEvent(input$commitCnvGroup, { commitCnvGroup() })

#----------------------------------------------------------------------
# handle clicks within chromosome-level views for marking CNV groups
#----------------------------------------------------------------------
initGenomePlotClicks <- reactiveVal(0)
observeEvent(initGenomePlotClicks(), {
    session$sendCustomMessage("cellPlotsWrapperInit", list(
        prefix = session$ns(""),
        divId = "chromPlotsWrapper"
    ))
    session$sendCustomMessage("cellPlotsWrapperUpdate", list(
        prefix = session$ns(""),
        divId = "chromPlotsWrapper",
        cellsPerPage = nWorkingCells(),
        short = TRUE
    ))
}, ignoreInit = TRUE)
observeEvent(input$cellWindowsPlotClick, {
    handleCellPlotClick(input$cellWindowsPlotClick, zoomChrom, zoomTargetWindow)
})
zoomTargetWindow <- reactiveVal(NULL)
pendingIndexKey <- NULL
observeEvent(zoomTargetWindow(), { # worker function called when a CNV area is clicked to (un)set its group
    click <- zoomTargetWindow()
    req(click)
    cnvs <- cnvsTableData()
    cnv <- cnvs[
        sourceId == click$click$data$source_id & 
        cell_id == as.character( click$click$data$cell_id ) &
        chrom == click$targetWindow$chrom & 
        click$targetWindow$start >= start & 
        click$targetWindow$start <= end
    ]
    if(nrow(cnv) != 1) return(NULL)

    assignments <- cnvGroupAssignments()$table
    cnvHasCommittedAssignment <- assignments[key == cnv$key & !is.na(groupKey), .N > 0] 
    if(cnvHasCommittedAssignment){ # if a CNV is already assigned, act on its current group in priority over any pending group
        if(assignments[key == cnv$key, key == groupKey]){ # delete an entire previously commited group
            inGroup <- assignments[, !is.na(groupKey) & groupKey == cnv$key]
            inGroupKeys <- assignments[inGroup, key]
            assignments <- assignments[!inGroup]
            invalidateZoomedCells <<- cnvs[key %in% inGroupKeys, paste(sourceId, cell_id)]
        } else { # remove a CNV from a committed group
            assignments <- assignments[key != cnv$key]
            invalidateZoomedCells <<- paste(click$click$data$source_id, click$click$data$cell_id)
        }
    } else {
        isPendingGroup <- !is.null(pendingIndexKey)
        if(isPendingGroup && cnv$key == pendingIndexKey){ # drop all pending assignments when the index of a pending group is clicked
            wasPending <- assignments[, !is.na(pendingKey)]
            wasPendingKeys <- assignments[wasPending, key]
            assignments[wasPending, pendingKey := as.character(NA)]
            assignments <- assignments[!is.na(groupKey)]
            invalidateZoomedCells <<- cnvs[key %in% wasPendingKeys, paste(sourceId, cell_id)] # update all affected cells
            pendingIndexKey <<- NULL
        } else {
            cnvHasAssignment <- assignments[key == cnv$key, .N > 0] 
            if(cnvHasAssignment && assignments[key == cnv$key, !is.na(pendingKey)]){ # remove a CNV from a pending group
                assignments[key == cnv$key, pendingKey := as.character(NA)]
                assignments <- assignments[!is.na(groupKey) | !is.na(pendingKey)]  
            } else {
                if(!isPendingGroup) pendingIndexKey <<- cnv$key # assign the index CNV for a new group as the first one clicked
                if(cnvHasAssignment){ # tentatively add a previously assigned CNV to the new group
                    assignments[key == cnv$key, pendingKey := pendingIndexKey]
                } else {
                    newRow <- data.table( # tentatively add a previously unassigned CNV to the new group
                        groupKey = as.character(NA),
                        pendingKey = pendingIndexKey
                    )
                    newRow[, key := cnv$key]
                    assignments <- rbind(assignments, newRow)
                }
            }
            invalidateZoomedCells <<- paste(click$click$data$source_id, click$click$data$cell_id) # update this one affected cell
        } 
    }
    cnvGroupAssignments(list(table = assignments, entropy = sample(1e8, 1)))
})
isMatchedCnv <- function(cnvKey){ # returns the display properties to plotCellQC()
    cnvs <- cnvTableWithGroups()
    cnv <- cnvs[key == cnvKey]  
    isPendingGroup <- !is.null(pendingIndexKey)
    hasPendingAssignment <- !is.na(cnv$pendingKey)
    isAssigned <- hasPendingAssignment || !is.na(cnv$groupKey)
    isIndexCnv <- cnv$key == (if(hasPendingAssignment) pendingIndexKey else cnv$groupKey)
    CN_col <- if(cnv$cellCN > cnv$referenceCN) 3 else 1
    list(
        isAssigned = isAssigned,
        CN_col = if(!isAssigned) 0 else CN_col,
        CN_border = if(!isAssigned) NA else if(hasPendingAssignment) 0 else CN_col,
        lty = if(!isAssigned) NULL else if(isIndexCnv) 1 else 2
    )
}
commitCnvGroup <- function(){ # commit a pending CNV group as a final outcome of this module
    req(!is.null(pendingIndexKey))
    assignments <- cnvGroupAssignments()$table
    wasPending <- assignments[, !is.na(pendingKey)]
    wasPendingKeys <- assignments[wasPending, key]
    assignments[wasPending, ":="(
        groupKey = pendingKey,
        pendingKey = as.character(NA)
    )]
    assignments <- assignments[!is.na(groupKey)]    
    cnvs <- cnvsTableData()
    invalidateZoomedCells <<- cnvs[key %in% wasPendingKeys, paste(sourceId, cell_id)] 
    pendingIndexKey <<- NULL
    cnvGroupAssignments(list(table = assignments, entropy = sample(1e8, 1)))
}

#----------------------------------------------------------------------
# handle clearing of group assignment under various triggers
#----------------------------------------------------------------------
observeEvent(zoomChrom(), { # clear pending assignments when switching chromosomes
    assignments <- cnvGroupAssignments()$table
    assignments[, pendingKey := as.character(NA)]
    assignments <- assignments[!is.na(groupKey)]
    pendingIndexKey <<- NULL
    cnvGroupAssignments(list(table = assignments, entropy = sample(1e8, 1))) 
})
unlinkProjectZoomPngs <- function(sourceIds, zoomChrom = "*"){
    for(sourceId in sourceIds){
        project <- getScCnvProjectData(sourceId)
        filename <- paste0("*.", zoomChrom, ".short.*.png")
        unlink(file.path(project$qcPlotsDir, filename))
    }
}
clearSourceAssignments <- function(zoomChrom = NULL){
    assignments <- cnvGroupAssignments()$table
    req(assignments)
    cnvs <- cnvsTableData()
    if(!is.null(zoomChrom)) cnvs <- cnvs[chrom == zoomChrom]
    assignments <- assignments[!(key %in% cnvs$key)]
    pendingIndexKey <<- NULL
    unlinkProjectZoomPngs(cnvs$sourceId, zoomChrom)
    cnvGroupAssignments(list(table = assignments, entropy = sample(1e8, 1)))       
}
observeEvent(input$clearChromAssignments, { # respond to special links to execute bulk clearing of assignments
    clearSourceAssignments(zoomChrom())
})
observeEvent(input$clearSourceAssignments, { # clear all chromosomes for the selected source(s)
    clearSourceAssignments()
})
observeEvent(input$clearAllAssignments, { # completely purge all assignment outcomes
    assignments <- cnvGroupAssignments()$table[FALSE]
    pendingIndexKey <<- NULL
    sources <- getStepReturnValueByType('upload', 'outcomes')$sources()
    unlinkProjectZoomPngs(names(sources))
    cnvGroupAssignments(list(table = assignments, entropy = sample(1e8, 1))) 
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    if(!is.null(bm$outcomes$assignments)) cnvGroupAssignments(bm$outcomes$assignments)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(
        assignments = cnvGroupAssignments
    ),
    addCnvGroupAssignments = function(cnvs){
        req(cnvs)
        assignments <- copy(cnvGroupAssignments()$table)
        if(nrow(assignments) > 0) {
            cnvs <- merge(
                cnvs, 
                assignments[, .SD, .SDcols = c("key", "groupKey")], 
                by = "key", all.x = TRUE, all.y = FALSE
            )
            cnvs[is.na(groupKey), groupKey := key]
        } else cnvs[, groupKey := key]
        cnvs
    },
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
