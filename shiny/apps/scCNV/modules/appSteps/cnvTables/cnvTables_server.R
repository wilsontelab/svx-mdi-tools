#----------------------------------------------------------------------
# server components for the cnvTables appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
cnvTablesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'cnvTables'
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
cnvsTableData <- reactive({ # a table of CNVs that were kept from the selected sources
    getFilteredCnvData(input, sourceIds, settings)
})
cnvTableWithGroups <- reactive({ # combine cnvs and group assignments into one working table
    cnvs <- cnvsTableData()
    req(cnvs)
    cnvs <- app$collateCNVs$addCnvGroupAssignments(cnvs)[order(chromI, start, end)]
    cnvs[, ":="(
        regions = paste0(chrom, ":", start, "-", end),
        cellKey = paste(sourceId, cell_id)
    )]

    # merge CNVs grouped within one cell into a single event row
    outputCnvs <- cnvs[order(chromI, start, end), .( 
        sourceId = sourceId[1],
        project = project[1],
        sampleName = sampleName[1],
        cell_id = cell_id[1],
        chrom = chrom[1],
        start = min(start, na.rm = TRUE),
        end = max(end, na.rm = TRUE),
        size = commify(max(end, na.rm = TRUE) - min(start, na.rm = TRUE) + 1),
        aneuploid = any(aneuploid),
        telomeric = any(telomeric),
        cellCN = cellCN[1],
        referenceCN = referenceCN[1],
        nWindows = sum(nWindows),
        windowPower = windowPower[1],
        regions = paste(regions, collapse = ","),
        groupKey = groupKey[1],
        chromI = chromI[1], 
        fused = .N > 1
    ), by = c("cellKey", "groupKey")]

    # record simplified information on the index CNV for each group
    setkey(cnvs, key)
    outputCnvs[, matchedTo := cnvs[groupKey, regions]]
    groupKeys <- unique(outputCnvs$groupKey)
    groupIs <- as.list(1:length(groupKeys))
    names(groupIs) <- groupKeys
    outputCnvs[, grouped := .N > 1, by = "groupKey"]
    outputCnvs[, group := unlist(groupIs[groupKey])]
    outputCnvs
})

#----------------------------------------------------------------------
# create a display table of all kept CNVs across all selected data sources, one row per fused CNV
#----------------------------------------------------------------------
displayCols <- c(
    "group",        
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
    "windowPower",
    "regions",
    "matchedTo",
    "fused",
    "grouped"
)
cnvsTableData_byCnv <- reactive({
    cnvs <- cnvTableWithGroups()[order(chromI, -aneuploid, start, end, cellCN - referenceCN), .SD, .SDcols = c(
        "sourceId",
        displayCols
    )]
})
cnvsTable <- bufferedTableServer(
    "cnvsTable",
    id,
    input,
    tableData = function() cnvsTableData_byCnv()[, .SD, .SDcols = displayCols],
    # editBoxes = list(),
    selection = 'none',
    # selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# create a display table of all kept CNVs across all selected data sources, one row per grouped CNV
#----------------------------------------------------------------------
cnvsTableData_byGroup <- reactive({
    cnvsTableData_byCnv()[, {
        i <- which.min(windowPower) # use the highest resolution cell as the reference (in case the user did not)
        samples <- getSampleLabels(input, sourceIds, project, sampleName)
        .(
            nSamples = length(unique(samples)),        
            nCells = .N,
            samples = paste(samples, collapse = " || "),
            regions = paste(regions, collapse = " || "),
            chrom = chrom[i],
            start = start[i],
            end = end[i],
            size = size[i],
            aneuploid = any(aneuploid),
            telomeric = any(telomeric),
            cellCN = cellCN[i],
            referenceCN = referenceCN[i],
            maxNWindows = max(nWindows),
            minWindowPower = min(windowPower),
            anyFused = any(fused)
        )
    }, by = "group"]
})
cnvsTable <- bufferedTableServer(
    "groupsTable",
    id,
    input,
    tableData = cnvsTableData_byGroup,
    # editBoxes = list(),
    selection = 'none',
    # selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# create a pivot table of all kept CNVs across all selected data sources, one row per grouped CNV, one column per sample, nCells in cells
#----------------------------------------------------------------------
nKeptCells <- reactive({
    app$keepReject$getNKeptCells(settings)   
})
cnvsTableData_pivot <- reactive({
    cnvsTableData_byCnv <- cnvsTableData_byCnv()
    req(cnvsTableData_byCnv)
    nKeptCells <- nKeptCells()
    groups <- cnvsTableData_byCnv[, {
        i <- which.min(windowPower) # use the highest resolution cell as the reference (in case the user did not)
        .(       
            chrom = chrom[i],
            start = start[i],
            end = end[i],
            size = size[i],
            aneuploid = any(aneuploid),
            telomeric = any(telomeric),
            cellCN = cellCN[i],
            referenceCN = referenceCN[i],
            sampleLabel = getSampleLabels(input, sourceIds, project, sampleName),
            nCells = .N,
            nKeptCells = nKeptCells[paste(sourceId, sampleName, sep = "::"), N]
        )
    }, by = c("group")]
    nCells <- dcast(
        groups, 
        chrom + start + end + size + aneuploid + telomeric + cellCN + referenceCN ~ sampleLabel, 
        value.var = "nCells",
        fun.aggregate = length
    )
    nKeptCells <- dcast(
        groups, 
        chrom + start + end + size + aneuploid + telomeric + cellCN + referenceCN ~ sampleLabel, 
        value.var = "nKeptCells",
        fun.aggregate = first
    )
    for(j in 9:ncol(nCells)) nCells[[j]] <- ifelse(nCells[[j]] == 0, "", paste(nCells[[j]], nKeptCells[[j]], sep = " / "))
    chromIs <- 1:100
    names(chromIs) <- paste0("chr", c(as.character(1:98), "X", "Y"))
    nCells[order(chromIs[chrom], start, end)]
})
cnvsTable <- bufferedTableServer(
    "groupsPivotTable",
    id,
    input,
    tableData = cnvsTableData_pivot,
    # editBoxes = list(),
    selection = 'none',
    # selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
