#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# help manage the server (not session) scCNV data cache
getScCnvProjectData <- function(sourceId){ # sourceId is a string key
    filePath <- loadPersistentFile(
        sourceId = sourceId, 
        contentFileType = "normalizeFile", 
        ttl = 10 * 60,
        postProcess = function(x){
            setkey(x$colData, "cell_id")  
            x$qcPlotsDir <- expandSourceFilePath(sourceId, "qc_plots")
            if(!dir.exists(x$qcPlotsDir)) dir.create(x$qcPlotsDir)
            x$rowRanges <- NULL # remove the objects that inflate the size of the full project file for faster caching
            x$raw_counts <- NULL
            x
        },
        session = session,
        spinnerMessage = "loading project data"
    )
    persistentCache[[filePath]]$data
}
scCnvProjectReactive <- function(sourceId, ...){
    reactive({
        if(is.function(sourceId)) sourceId <- sourceId()
        req(sourceId)
        getScCnvProjectData(sourceId)
    }) 
}
getScCnvCellMatrix <- function(sourceId, settings, minCellWindowPower = NULL, col = "CN"){
    project <- getScCnvProjectData(sourceId)
    if(is.null(minCellWindowPower)) minCellWindowPower <- project$colData[bad == FALSE, min(windowPower)]
    shapeModel <- getShapeModel(settings, cells = project$cells) 
    repKey <- "sequential"
    matrixFile <- paste("cell_matrix", shapeModel$key, repKey, minCellWindowPower, "rds", sep = ".")
    filePath <- loadPersistentFile(
        file = expandSourceFilePath(sourceId, matrixFile),
        force = FALSE,
        ttl = 60 * 60,
        create = function(file){
            I <- !project$colData$bad
            cell_ids <- project$colData[I, cell_id]
            m <- sapply(cell_ids, function(cell_id){
                cell <- project$cells[[cell_id]]
                x <- cell$windows[[shapeModel$key]][[repKey]][[col]]
                if(cell$windowPower > minCellWindowPower){
                    expandBy <- 2 ** (cell$windowPower - minCellWindowPower)
                    x <- project$windows[[cell$windowPower + 1]][, {
                        chrom_ <- chrom
                        nExpandedWindows <- project$windows[[minCellWindowPower + 1]][chrom == chrom_, .N]
                        uncollapseVector(x[.I], expandBy, nExpandedWindows)
                    }, by = chrom][[2]]
                }
                x
            })
            colnames(m) <- cell_ids
            saveRDS(m, file = file)            
        },
        session = session,
        spinnerMessage = "loading cell matrix"
    )
    persistentCache[[filePath]]$data
}
getAllSourceCnvs <- function(sourceId, settings){ # get all CNVs called by the pipeline for a specific data source
    project <- getScCnvProjectData(sourceId)
    shapeModel <- getShapeModel(settings, cells = project$cells) 
    cnvsFile <- paste("source_cnvs", shapeModel$key, "rds", sep = ".")
    filePath <- loadPersistentFile(
        file = expandSourceFilePath(sourceId, cnvsFile),
        force = FALSE,
        ttl = 60 * 60,
        create = function(file){
            cnvs <- project$cnvs[[shapeModel$key]]
            N <- if(is.null(cnvs)) 0 else nrow(cnvs)
            x <- if(N == 0) NULL else {
                cnvs[, key := getCnvKey(.SD, sourceId), by = seq_len(N)]
                cnvs
            }
            saveRDS(x, file = file)            
        },
        session = session,
        spinnerMessage = "loading source cnvs"
    )
    persistentCache[[filePath]]$data
}

# help manage the session cache (not server, depends on user actions)
getCnvsTableData <- function(sourceIds, settings){ # a table of CNVs that were kept from the selected sources
    req(sourceIds)
    keptCnvKeys <- app$keepReject$getKeptCnvKeys()
    req(keptCnvKeys)
    chromIs <- 1:100
    names(chromIs) <- paste0("chr", c(as.character(1:98), "X", "Y"))

    # extend and concatenate CNVs from all selected sources
    do.call(rbind, lapply(sourceIds, function(sourceId){
        cnvs <- getAllSourceCnvs(sourceId, settings)
        if(is.null(cnvs)) return(NULL)
        cnvs <- cnvs[key %in% keptCnvKeys] 
        project <- getScCnvProjectData(sourceId)
        chromRanges <- project$windows[[length(project$windows)]][, .(
            start = min(start, na.rm = TRUE),
            end   = max(end, na.rm = TRUE)
        ), by = "chrom"]
        getWStartI <- Vectorize(function(cell_id, cnvChrom, cnvStart){
            w <- project$windows[[project$cells[[cell_id]]$windowPower + 1]]
            w[chrom == cnvChrom & start ==cnvStart, i]
        })
        cnvs[, ":="(
            sourceId = sourceId,
            project = getSourceFilePackageName(sourceId),
            chromI = chromIs[chrom],
            aneuploid = mapply(function(chrom_, start_, end_){
                chromRanges[chrom == chrom_, start == start_ && end == end_]
            }, chrom, start, end),
            telomeric = mapply(function(chrom_, start_, end_){
                chromRanges[chrom == chrom_, start == start_ || end == end_]
            }, chrom, start, end),
            size = commify(end - start + 1),
            CN = ifelse(cellCN > referenceCN, 3, 1),
            startI = getWStartI(cell_id, chrom, start)
        )]
        cnvs[, endI := startI + nWindows - 1]
        cnvs
    }))[order(chromI, -aneuploid, start, end, cellCN - referenceCN)] # by chromosome, aneuploid first in chrom group
}

# shared support for certain appStep modules UI
ALL <- "__ALL__"
scCnvDataSourceTableUI <- function(ns){
    dataSourceTableUI(
        ns("source"), 
        "Project Source", 
        width = 8, 
        collapsible = TRUE
    ) 
}
sampleSelectorDivUI <- function(ns){
    tags$div(
        style = "margin-top: 0; margin-bottom: 10px; white-space: nowrap;",
        tags$div(
            class = "cellPageInput",
            style = "width: 400px; margin-right: 10px;", 
            selectInput(ns('sampleNameFilter'), "Sample", choices = c(All = "__ALL__"), width = "100%")
        )
    )
}
updateSampleSelector <- function(session, input, sourceIds){
    x <- unique(do.call(rbind, lapply(sourceIds, function(sourceId){
        project <- getScCnvProjectData(sourceId)
        data.table(
            sourceId = sourceId, 
            Project = project$manifest$Project,
            Sample_Name = project$manifest$Sample_Name
        )
    })))
    choices <- x[, paste(sourceId, Sample_Name, sep = "::")]
    names(choices) <- x[, if(length(sourceIds) == 1) Sample_Name else paste(Project, Sample_Name, sep = " / ")]
    freezeReactiveValue(input, "sampleNameFilter")
    updateSelectInput(session, "sampleNameFilter", choices = c(All = ALL, choices), selected = ALL)    
}
getFilteredCnvData <- function(input, sourceIds, settings){
    req(input$sampleNameFilter)
    cnvs <- getCnvsTableData(sourceIds(), settings)
    cnvs[if(input$sampleNameFilter == ALL) TRUE else paste(sourceId, sampleName, sep = "::") == input$sampleNameFilter]
}
getSampleLabels <- function(input, sourceIds, projects, sampleNames){
    isSingleSource <- length(sourceIds()) == 1 || input$sampleNameFilter != ALL  
    if(isSingleSource) sampleNames else ifelse(
        projects == sampleNames,
        sampleNames,
        paste(projects, sampleNames, sep = "/") 
    ) 
}

# appServer function called after all modules are instantiated
appServer <- function(){

}
