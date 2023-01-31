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
matchedCnvPairs <- reactiveValues()

#----------------------------------------------------------------------
# load sample/cell source data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "multiple"
)
allSourceCnvs <- list()
getAllSourceCnvs <- function(sourceId){
    if(is.null(allSourceCnvs[[sourceId]])){
        startSpinner(session, message = "loading source CNVs")
        project <- getScCnvProjectData(sourceId)
        shapeModel <- getShapeModel(settings, cells = project$cells)        
        cnvs <- project$cnvs[[shapeModel$key]]
        N <- if(is.null(cnvs)) 0 else nrow(cnvs)
        allSourceCnvs[[sourceId]] <<- if(N == 0) NULL else {
            cnvs[, key := getCnvKey(.SD, sourceId), by = seq_len(N)]
            cnvs
        }
        stopSpinner(session)
    }
    allSourceCnvs[[sourceId]]
}
cnvsTableData <- reactive({
    sourceIds <- sourceId()
    req(sourceIds)
    keptCnvKeys <- app$keepReject$getKeptCnvKeys()
    req(keptCnvKeys)
    chromIs <- 1:100
    names(chromIs) <- paste0("chr", c(as.character(1:98), "X", "Y"))

    # extend and concatenate CNVs from all selected sources
    cnvs <- do.call(rbind, lapply(sourceIds, function(sourceId){
        cnvs <- getAllSourceCnvs(sourceId)
        if(is.null(cnvs)) return(NULL)
        cnvs <- cnvs[key %in% keptCnvKeys] 
        project <- getScCnvProjectData(sourceId)
        chromRanges <- project$windows[[length(project$windows)]][, .(
            start = min(start, na.rm = TRUE),
            end   = max(end, na.rm = TRUE)
        ), by = "chrom"]
        cnvs[, ":="(
            project = getSourceFilePackageName(sourceId),
            chromI = chromIs[chrom],
            aneuploid = mapply(function(chrom_, start_, end_){
                chromRanges[chrom == chrom_, start == start_ && end == end_]
            }, chrom, start, end),
            telomeric = mapply(function(chrom_, start_, end_){
                chromRanges[chrom == chrom_, start == start_ || end == end_]
            }, chrom, start, end),
            size = commify(end - start + 1)
        )]
        cnvs
    }))[order(chromI, -aneuploid, start, end, cellCN - referenceCN)] # by chromosome, aneuploid first in chrom group

    # find overlap groups of CNVs
    # TODO: should overlap group use CNC?
    # TODO: manual matching...
    overlapGroup <- 0
    incrementOverlapGroup <- function(){
        overlapGroup <<- overlapGroup + 1
        overlapGroup
    }
    getChromOverlapGroups <- function(cnvs){
        if(nrow(cnvs) == 1) return(incrementOverlapGroup())
        nAneuploid <- cnvs[, sum(aneuploid)]
        overlapGroups <- if(nAneuploid > 0) rep(incrementOverlapGroup(), nAneuploid) else integer()
        if(cnvs[, all(aneuploid)]) return(overlapGroups)
        i <- nAneuploid + 1
        overlapGroups <- c(overlapGroups, incrementOverlapGroup())
        if(i == nrow(cnvs)) return(overlapGroups)
        for(j in (i + 1):nrow(cnvs)){
            if(cnvs[j - 1, start] <= cnvs[j,     end] && 
               cnvs[j,     start] <= cnvs[j - 1, end]){
                overlapGroups <- c(overlapGroups, overlapGroups[length(overlapGroups)])
            } else {
                overlapGroups <- c(overlapGroups, incrementOverlapGroup())
            }  
        }
        overlapGroups        
    }
    cnvs[, overlapGroup := getChromOverlapGroups(.SD), by = "chrom"]
    cnvs
})
cnvsTableDataPretty <- reactive({
    cnvs <- cnvsTableData()[, .SD, .SDcols = c(
        "overlapGroup",
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

#----------------------------------------------------------------------
# create a table of all kept CNVs across all selected data sources
#----------------------------------------------------------------------
cnvsTable <- bufferedTableServer(
    "cnvsTable",
    id,
    input,
    tableData = cnvsTableDataPretty,
    # editBoxes = list(),
    selection = 'single',
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
