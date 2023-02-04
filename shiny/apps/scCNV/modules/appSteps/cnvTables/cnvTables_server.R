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
# load sample/cell source data
#----------------------------------------------------------------------
sourceIds <- dataSourceTableServer(
    "source", 
    selection = "multiple"
)
cnvsTableData <- reactive({ # a table of CNVs that were kept from the selected sources
    getCnvsTableData(sourceIds(), settings)
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
cnvsTableData_byCnv <- reactive({
    cnvs <- cnvTableWithGroups()[order(chromI, -aneuploid, start, end, cellCN - referenceCN), .SD, .SDcols = c(
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
    )]
})
cnvsTable <- bufferedTableServer(
    "cnvsTable",
    id,
    input,
    tableData = cnvsTableData_byCnv,
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
        samples <- paste(project, sampleName)
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
cnvsTableData_pivot <- reactive({
    groups <- cnvsTableData_byCnv()[, {
        i <- which.min(windowPower) # use the highest resolution cell as the reference (in case the user did not)
        .(       
            nCells = .N,
            sample = paste(project, sampleName, sep = "/"),
            chrom = chrom[i],
            start = start[i],
            end = end[i],
            size = size[i],
            aneuploid = any(aneuploid),
            telomeric = any(telomeric),
            cellCN = cellCN[i],
            referenceCN = referenceCN[i]
        )
    }, by = c("group")]
    dcast(
        groups, 
        group + chrom + start + end + aneuploid + telomeric + cellCN + referenceCN ~ sample, 
        value.var = "nCells"
    )
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
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
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
