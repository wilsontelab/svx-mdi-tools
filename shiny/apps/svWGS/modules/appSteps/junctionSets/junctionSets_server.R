#----------------------------------------------------------------------
# server components for the junctionSets appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
junctionSetsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'junctionSets'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
gridColors <- list(x = "#aaaaaa", y = "#aaaaaa")

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single") 
sampleNames <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    source <- getSourceFromId(sourceId)
    req(source)
    source$manifest$Sample_ID
})
genomes <- reactive({ # data.table with genome,binSize,isComposite,compositeDelimiter
    sourceId <- sourceId()
    req(sourceId)
    readRDS(getSourceFilePath(sourceId, "genomesFile"))
})
junctions <- reactive({
    sourceId <- sourceId()
    req(sourceId)  
    x <- svx_loadJunctions(sourceId, svWGS_loadJunctions) 
    stopSpinner(session)
    x
})
getGenomeJxns <- function(genome, junctions = NULL) {
    if(is.null(junctions)) junctions <- junctions()
    if(genome$isComposite) junctions[
        endsWith(CHROM_1, genome$genome) & # allow translocations but suppress intergenomic
        endsWith(CHROM_2, genome$genome)
    ] else junctions
}

#----------------------------------------------------------------------
# source samples, stratified by composite genomes when applicable
#----------------------------------------------------------------------
sampleGenomeTableData <- reactive({
    as.data.table(expand.grid(
        Sample = sampleNames(), 
        Genome = genomes()$genome, 
        stringsAsFactors = FALSE
    ))[, ":="(
        key = paste(Sample, Genome),
        sourceId = sourceId()        
    )]
})
sampleGenomeTable <- bufferedTableServer(
    "sampleGenome",
    id,
    input,
    tableData = reactive( sampleGenomeTableData()[, .(Sample, Genome)] ),
    selection = 'single',
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
)
sampleGenome <- reactive({
    row <- sampleGenomeTable$rows_selected()
    req(row)
    sampleGenomeTableData()[row]
})

#----------------------------------------------------------------------
# gc bias models read from NormalizeGC
#----------------------------------------------------------------------
gcBiasModels <- reactive({
    app$normalizeGC$getGcBiasModels(sourceId())
})
gcBiasModel <- reactive({
    gcBiasModels <- gcBiasModels()
    sampleGenome <- sampleGenome()
    gcBiasModels[[sampleGenome$key]]
})

#----------------------------------------------------------------------
# interactive plot to correlate two junction properties
#----------------------------------------------------------------------
cnAxisLim  <- c(0,4)
cncAxisLim <- c(-3.5,3.5)
junctionAxisChoices <- list(
    junctionCN = list(
        lab = "Junction CN",
        lim = cnAxisLim
    ),
    outerFlankCN = list(
        lab = "Outer Flank CN",
        lim = cnAxisLim
    ),
    innerFlankCN = list(
        lab = "Inner Flank CN",
        lim = cnAxisLim
    ),
    eventSpanCN = list(
        lab = "Event Span CN",
        lim = cnAxisLim
    ),
    flankCNC = list(
        lab = "Flank CNC",
        lim = cncAxisLim
    ),
    eventCNC = list(
        lab = "Event Span CNC",
        lim = cncAxisLim
    ),
    chromCNC = list(
        lab = "CNC",
        lim = cncAxisLim
    ),
    N_SAMPLES = list(
        lab = "# of Samples",
        lim = c(0.5, 2.5)
    ),
    SV_SIZE = list(
        lab = "log10(SV Size)",
        lim = c(2.9, 7.5)
    )
)
jxnValue <- function(d, col) switch(
    col,
    SV_SIZE = log10(d[[col]]),
    N_SAMPLES = jitter(d[[col]]),
    d[[col]]
)
jxnPlotData <- reactive({
    gcBiasModel <- gcBiasModel()
    if(!isTruthy(gcBiasModel)) return(NULL) 
    edgeTypes <- svx_jxnType_nameToX(input$jxnPlotSVTypes, "code")
    d <- gcBiasModel$jxnCN[edgeType %in% edgeTypes &
                          N_TOTAL  >= input$jxnPlotMinNInstances &
                          N_SPLITS >= input$jxnPlotMinSequenced]
    d <- if(input$jxnPlotNSamples == "Unique") d[N_SAMPLES == 1] # therefore, restricts to junctions unique to this sample
         else d[N_SAMPLES > 1]                                   # therefore, excludes     junctions unique to this sample
    d[, ":="(
        x = jxnValue(d, input$jxnPlotXAxis),
        y = jxnValue(d, input$jxnPlotYAxis),
        color = svx_jxnType_codeToX(d$edgeType, "color")
    )]
    setkey(d, svId)
    d
})
jxnPlot <- interactiveScatterplotServer(
    "jxnPlot",
    plotData = jxnPlotData,
    accelerate = TRUE,
    pointSize = 4,
    xtitle = reactive( junctionAxisChoices[[input$jxnPlotXAxis]]$lab),
    xrange = reactive( junctionAxisChoices[[input$jxnPlotXAxis]]$lim),
    ytitle = reactive( junctionAxisChoices[[input$jxnPlotYAxis]]$lab),
    yrange = reactive( junctionAxisChoices[[input$jxnPlotYAxis]]$lim),
    grid = gridColors,
    selectable = TRUE,
    keyColumn = "svId"
)  
jxnPlotSelected <- reactive({
    selected <- jxnPlot$selected() 
    req(selected, nrow(selected) > 0)    
    jxnPlotData()[selected$customdata] # as keyed by svId, passed to custom data via keyColumn above    
})

#----------------------------------------------------------------------
# interactive plot to correlate two junction properties, applied to the gated susbset of junctions from plot above
#----------------------------------------------------------------------
gatedJxnPlotData <- reactive({
    d <- jxnPlotSelected()
    d[, ":="(
        x = jxnValue(d, input$gatedJxnPlotXAxis),
        y = jxnValue(d, input$gatedJxnPlotYAxis)
    )]
    setkey(d, svId)
    d
})
gatedJxnPlot <- interactiveScatterplotServer(
    "gatedJxnPlot",
    plotData = gatedJxnPlotData,
    accelerate = TRUE,
    pointSize = 4,
    xtitle = reactive( junctionAxisChoices[[input$gatedJxnPlotXAxis]]$lab),
    xrange = reactive( junctionAxisChoices[[input$gatedJxnPlotXAxis]]$lim),
    ytitle = reactive( junctionAxisChoices[[input$gatedJxnPlotYAxis]]$lab),
    yrange = reactive( junctionAxisChoices[[input$gatedJxnPlotYAxis]]$lim),
    grid = gridColors,
    selectable = TRUE,
    keyColumn = "svId"
)  
gatedJxnPlotSelected <- reactive({
    selected <- gatedJxnPlot$selected() 
    if(!isTruthy(selected)) return(jxnPlotSelected())
    gatedJxnPlotData()[selected$customdata] # as keyed by svId, passed to custom data via keyColumn above    
})

#----------------------------------------------------------------------
# table of filtered and gated SVs
#----------------------------------------------------------------------
gatedSvsData <- reactive({
    gcBiasModel <- gcBiasModel()
    req(gcBiasModel, gcBiasModel$genome)    
    jxns <- getGenomeJxns(gcBiasModel$genome)
    selected <- gatedJxnPlotSelected()
    jxns[selected$svId]
})
gatedSvs <- bufferedTableServer(
    "gatedSvs",
    id,
    input,
    tableData = reactive({
        d <- gatedSvsData()
        d[, .(
            svId = SV_ID,
            type = svx_jxnType_codeToX(edgeType, "name"),
            size = SV_SIZE,
            insertSize,
            nInstances,
            nSequenced,
            flankLength,
            nLinkedJxns = nLinkedJunctions,
            rStart = paste0(cChrom1, ":", cRefPos1, ifelse(cStrand1 == 1, "+", "-")),
            rEnd   = paste0(cChrom2, ":", cRefPos2, ifelse(cStrand2 == 1, "+", "-"))
        )] 
    }),
    selection = 'single',
    options = list(
    )
)
gatedSvsSelected <- reactive({
    selected <- gatedSvs$rows_selected() 
    req(selected)
    gatedSvsData()[selected]
})

#----------------------------------------------------------------------
# expanded molecule support for an SV junction selected in gatedSvs table
#----------------------------------------------------------------------
observeEvent(gatedSvsSelected(), {
    svWGS_expandJunction_(
        jxn = gatedSvsSelected(),
        targetId = sourceId(),
        objectTableFn    = NULL, 
        expansionTableFn = expansionTableData, 
        expansionUIFn    = expansionUIContents   
    )
    stopSpinner(session)

    # jxn <- 

    # # TODO: used the new svWGS_expandJunction_

    # molecules <- svWGS_loadMolecules(sourceId(), jxn$SV_ID)
    # molecules %>% 
    # svWGS_expansionTable(jxn) %>% 
    # expansionTableData()

    # # if the junction was sequenced, create the expansion2 elements
    # # a map and an alignment at base-level detail
    # junctionMap <- tryCatch({
    #     getJunctionMap(list(sv = jxn, mols = molecules[sample.int(.N)]))
    # }, error = function(e) {
    #     stopSpinner(session)
    #     NULL
    # })
    # if(is.null(junctionMap)){
    #     expansionUIContents("")
    # } else {
    #     startSpinner(session, message = "analyzing junction")
    #     expansionUIContents(tagList(
    #         tryCatch({
    #             junctionMapTrackExpansionUI(NULL, junctionMap, Pixels_Per_Base = 2) 
    #         }, error = function(e) ""),
    #         tryCatch({ 
    #             junctionAlignmentTrackExpansionUI(NULL, junctionMap, Bases_Per_Line = 100, Alignment_Mode = "Reference Molecule") 
    #         }, error = function(e) { print(e); "" })
    #     ))     
    # }
})
expansionTableData <- reactiveVal(NULL)
expansionTable <- bufferedTableServer(
    "expansionTable",
    id,
    input,
    tableData = expansionTableData,
    selection = 'none',
    options = list(
    )
)
expansionUIContents <- reactiveVal()
output$expansionUI <- renderUI({
    expansionUIContents()
})

#----------------------------------------------------------------------
# saving gated SV sets
#----------------------------------------------------------------------
workingId <- NULL # set to a plot id when editing a previously saved set
sendFeedback <- function(x, ...) output$saveJunctionSetFeedback <- renderText(x)
getJunctionSetName <- function(id){
    name <- savedJunctionSets$names[[id]] # user name overrides
    if(is.null(name)) savedJunctionSets$list[[id]]$Name else name
}
getJunctionSetNames <- function(rows = TRUE){
    sapply(names(savedJunctionSets$list)[rows], getJunctionSetName)
}
savedJunctionSetsTemplate <- data.table(
    Remove          = character(),
    Name            = character(),
    Sample          = character(),
    Genome          = character(),
    SV_Types        = character(),
    Min_Instances   = character(),
    N_Samples       = character(),
    N_Junctions     = integer(),
    Hash            = character()
)
savedJunctionSets <- summaryTableServer(
    id = 'savedJunctionSets', # NOT ns(id) when nesting modules!
    parentId = id,
    stepNumber = options$stepNumber,
    stepLocks = locks[[id]],
    sendFeedback = sendFeedback,
    template = savedJunctionSetsTemplate,
    type = 'shortList',
    remove = list(
        message = "Remove this set of saved junctions?",
        name = getJunctionSetName
    ),
    names = list(
        get = getJunctionSetNames,
        source = id
    )
) 
observeEvent(input$saveJunctionSet, {
    sourceId <- sourceId()
    gcBiasModel <- gcBiasModel()
    jxns <- gatedJxnPlotSelected()
    req(sourceId, gcBiasModel, input$saveJunctionSet, input$saveJunctionSet != 0, jxns, nrow(jxns) > 0)
    d <- list( # plot-defining metadata, shown on Saved Plots table; these define the data available to the plot
        Name = paste("SV Set #", length(savedJunctionSets$list) + 1),
        # Source = getSourceFilePackageName(sourceId), # use the source name, not its unique ID, to allow sample additions to saved plots
        Sample          = gcBiasModel$sampleGenome$Sample,
        Genome          = gcBiasModel$sampleGenome$Genome,
        SV_Types        = input$jxnPlotSVTypes,
        Min_Instances   = input$jxnPlotMinNInstances,
        Min_Sequence    = input$jxnPlotMinSequenced,
        N_Samples       = input$jxnPlotNSamples,
        N_Junctions     = nrow(jxns),
        Hash            = digest(jxns)
    )
    r <- initializeRecordEdit(d, workingId, savedJunctionSets$list, 'Junction Set', 'junction set', sendFeedback)
    d <- c(d, list( # non-definining attributes saved with junction set but not displayed on Saved Plots table
        sourceId = sourceId(),
        junctions = jxns
    ))
    saveEditedRecord(d, workingId, savedJunctionSets, r)
    workingId <<- NULL
})
addJunctionSetColumns <- function(dt, r, name = FALSE){
    for(x in c("Sample","Genome","Min_Instances","N_Samples","N_Junctions","Hash")) dt[[x]] <- r[[x]]
    for(x in c("SV_Types")) dt[[x]] <- paste(r[[x]], collapse = " ")
    dt
}
addDataListObserver(module, savedJunctionSetsTemplate, savedJunctionSets, function(r, ...){
    dt <- data.table(
        Remove = '', 
        Name   = ''
    )
    addJunctionSetColumns(dt, r)
})

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    savedJunctionSets$list  <- bm$outcomes$junctionSets
    savedJunctionSets$names <- bm$outcomes$junctionSetNames
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(
        junctionSets     = reactive(savedJunctionSets$list),
        junctionSetNames = reactive(savedJunctionSets$names)
    ),
    getSavedJunctionSets = reactive({ # lists saved junction sets, used to populate the junctionSets track items list, etc.
        jss   <- app$junctionSets$outcomes$junctionSets()
        names <- app$junctionSets$outcomes$junctionSetNames()
        do.call(rbind, lapply(names(jss), function(id){
            r <- jss[[id]]            
            dt <- data.table(Name = if(is.null(names[[id]])) r$Name else names[[id]]) # user name overrides
            addJunctionSetColumns(dt, r)
        }))        
    }),
    getJunctionSetSvs = function(hashes){ # expand the junctions found in a sample set, for navigator tables, etc.
        jss <- app$junctionSets$outcomes$junctionSets()
        req(hashes, jss)   
        do.call(rbind, lapply(jss, function(js){
            if(js$Hash %in% hashes) {
                merge(
                    svx_loadJunctions(js$sourceId, svWGS_loadJunctions)[SV_ID %in% js$junctions$svId],
                    app$normalizeGC$getCnvJxnsNormalizedCN(js$sourceId)$value[, .SD, .SDcols = c("SV_ID", "outerFlankCN", "innerFlankCN", "flankCNC", "junctionCN")],
                    by = "SV_ID",
                    all.x = TRUE
                )
            } else NULL
        }))
    },
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
