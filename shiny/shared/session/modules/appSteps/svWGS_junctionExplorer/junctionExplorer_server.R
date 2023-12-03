#----------------------------------------------------------------------
# server components for the junctionExplorer appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svWGS_junctionExplorerServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'junctionExplorer'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    templates = list(file.path(app$sources$suiteSharedTypesDir, "browserTrackTypes/svFilter_settings.yml")), # a path to settings.yml, sent to settingsServer()
    # immediate = TRUE # plus any other arguments passed to
)
gridColors <- list(x = "#aaaaaa", y = "#aaaaaa")

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
data <- svWGS_explorer_dataSelectorsServer(
    "dataSelectors",
    id,
    settings, 
    sourceSelection = "multiple",
    sampleGenomesSelection = "multiple"
)

#----------------------------------------------------------------------
# microhomology/insertion distribution plot
#----------------------------------------------------------------------
parseGroupingCols <- function(jxns, lhsCols, groupingCols){
    droppedCols <- character()
    if(length(groupingCols) > 0) for(col in groupingCols){
        if(length(unique(jxns[[col]])) == 1) droppedCols <- c(droppedCols, col)
    }
    groupingCols <- groupingCols[!(groupingCols %in% droppedCols)]
    if(length(groupingCols) == 0) groupingCols <- "X"
    jxns[, group := .SD[, apply(.SD, 1, paste, collapse = ", "), .SDcols = groupingCols]]
    groups <- sort(unique(jxns$group))
    nGroupingCols <- length(groupingCols)
    groupCounts <- jxns[, .N, by = .(group)]
    setkey(groupCounts, group)
    list(
        groupingCols = groupingCols,
        nGroupingCols = length(groupingCols),
        groups = groups,
        nGroups = length(groups),
        hasGroups = nGroupingCols > 1 || groupingCols != "X",
        groupCounts = groupCounts[groups, N],
        jxns = jxns[, .SD, .SDcols = c(lhsCols, "SV_SIZE", "group")]
    )
}
dcastSvsByGroup <- function(svs, lhsCols, groups_, step = 1){
    dt <- dcast(
        svs,
        as.formula(paste(lhsCols, "~ group")), 
        value.var = names(svs)[1], 
        fun.aggregate = length, 
        fill = 0
    )
    for(gl in groups_) dt[[gl]] <- dt[[gl]] / sum(dt[[gl]], na.rm = TRUE)
    allBins <- seq(min(dt[[1]]), max(dt[[1]]), step)
    missingBins <- allBins[!(allBins %in% dt[[1]])]
    if(length(missingBins) > 0){
        missingBins <- data.table(tmp = missingBins)
        setnames(missingBins, "tmp", names(dt)[1])
        dt <- rbind(dt, missingBins, fill = TRUE)
        dt[is.na(dt)] <- 0
        setorderv(dt, names(dt)[1], order = 1L)
    } 
    dt[, .SD, .SDcols = c(lhsCols, groups_)]
}
microhomologyData <- reactive({
    jxns <- data$sampleGenomeJunctions()[N_SPLITS > 0]
    req(jxns, nrow(jxns) > 0)
    jxns[, insertSize := -MICROHOM_LEN]
    # lhsCols <- "insertSize"
    # grouping <- parseGroupingCols(
    #     jxns, 
    #     lhsCols,
    #     input$microhomologyGroupBy
    # )
    list(
        # grouping = grouping,
        # dists = dcastSvsByGroup(
        #     grouping$jxns, 
        #     lhsCols,
        #     grouping$groups, 
        #     step = 1
        # ),
        jxns = jxns
    )
})
#----------------------------------------------------------------------
insertSizeDensityPlot <- mdiDensityPlotBoxServer(
    id = "insertSizeDensityPlot",
    data = reactive( microhomologyData()$jxns ),
    lhsCols = "insertSize", # columns on the left hand side of the dcast formula, usually the X-axis value
    groupingCols = reactive( input$microhomologyGroupBy ), # columns used to define the groups to summarize
    step = 1, # the bin resolution on the X axis
    xlab = "Insert Size (bp)",
    extraCols = NULL # additional column to maintain from the original data
)

# microhomologyPlot <- staticPlotBoxServer(
#     "microhomologyPlot",
#     margins = TRUE,
#     title = TRUE,
#     settings = list(
#         Limits = list(
#             Min_X_Value = list(
#                 type = "numericInput",
#                 value = -10,
#                 min = -50, 
#                 max = 0,
#                 step = 5
#             ),
#             Max_X_Value = list(
#                 type = "numericInput",
#                 value = 15,
#                 min = 0, 
#                 max = 50,
#                 step = 5
#             )
#         )
#     ), 
#     size = "m",
#     create = function() {
#         d <- microhomologyData()
#         xlim <- c(
#             microhomologyPlot$settings$get("Limits","Min_X_Value"),
#             microhomologyPlot$settings$get("Limits","Max_X_Value")
#         )
#         maxY <- d$dists[, max(.SD, na.rm = TRUE), .SDcols = d$grouping$groups] * 1.05 
#         title <- microhomologyPlot$settings$get("Plot_Frame", "Title", NULL)
#         totalJxns <- paste(trimws(commify(sum(d$grouping$groupCounts))), "Total Junctions")
#         title <- if(is.null(title)) totalJxns else paste0(title, " (", totalJxns, ")")
#         microhomologyPlot$initializeFrame(
#             xlim = xlim,
#             ylim = c(0, maxY),
#             xlab = "Insert Size (bp)",
#             ylab = "Frequency",
#             xaxs = "i",
#             yaxs = "i",
#             title = title,
#             cex.main = 0.95
#         )
#         abline(v = c(seq(-50, 50, 5), -1, -2), col = "grey")
#         abline(v = 0) 
#         lwd <- 1.5
#         for(i in 2:ncol(d$dists)){ # overplot individual data points on the bar plot
#             microhomologyPlot$addLines(
#                 x = d$dists[[1]],
#                 y = d$dists[[i]],
#                 col = CONSTANTS$plotlyColors[[i - 1]],
#                 lwd = lwd
#             )
#         }
#         if(d$grouping$hasGroups) microhomologyPlot$addMarginLegend(
#             xlim[2] * 1.1, maxY, lty = 1, lwd = 1.5, 
#             legend = paste0(d$grouping$groups, " (", trimws(commify(d$grouping$groupCounts)), ")"),
#             col = unlist(CONSTANTS$plotlyColors[1:d$grouping$nGroups]),
#             bty = "n"
#         )
#         stopSpinner(session)
#     }
# )

#----------------------------------------------------------------------
# insert size vs. SV size
#----------------------------------------------------------------------
sizeCorrelationPlot <- svx_sizeCorrelationPlotBoxServer(
    "sizeCorrelationPlot",
    reactive( microhomologyData()$jxns )
)

#----------------------------------------------------------------------
# matching junctions table
#----------------------------------------------------------------------
matchingJunctions <- bufferedTableServer(
    "matchingJunctions",
    id,
    input,
    tableData = reactive({
        microhomologyData()$jxns[, .(
            svId = SV_ID,
            type = edgeType, 
            size = SV_SIZE,
            nTot = nInstances, 
            nSeq = nSequenced,
            nStr1 = STRAND_COUNT1, 
            nStr2 = STRAND_COUNT2,
            sample = SAMPLES,
            node1 = paste0(cChrom1, ":", cRefPos1, "/", cStrand1),
            node2 = paste0(cChrom2, ":", cRefPos2, "/", cStrand2),
            merge = MERGE_LEN, 
            jxnNBases = insertSize, 
            jxnBases =JXN_BASES
            # , 
            # jxnSeq = JXN_SEQ
        )] 
    }),
    selection = 'single',
    options = list()
)
matchingJunctionsSelected <- reactive({
    selected <- matchingJunctions$rows_selected() 
    req(selected)
    microhomologyData()$jxns[selected]
})
# [1] "SV_ID"                 "MAPQ_1"                "UMI_1"                
# mdi-web-server-app-server-1              |  [4] "N_TOTAL"               "N_GAPS"                "N_SPLITS"             
# mdi-web-server-app-server-1              |  [7] "N_OUTER_CLIPS"         "JXN_TYPE"              "N_DUPLEX"             
# mdi-web-server-app-server-1              | [10] "N_DUPLEX_GS"           "STRAND_COUNT"          "STRAND_COUNT_GS"      
# mdi-web-server-app-server-1              | [13] "STRAND_COUNT1"         "STRAND_COUNT2"         "TARGET_CLASS"         
# mdi-web-server-app-server-1              | [16] "SHARED_PROPER"         "SHARED_PROPER_GS"      "SAMPLES"              
# mdi-web-server-app-server-1              | [19] "N_SAMPLES"             "MAPQ_2"                "UMI_2"                
# mdi-web-server-app-server-1              | [22] "CHROM_1"               "SIDE_1"                "POS_1"                
# mdi-web-server-app-server-1              | [25] "CHROM_2"               "SIDE_2"                "POS_2"                
# mdi-web-server-app-server-1              | [28] "JUNCTION_NAME"         "JUNCTION_NAMES"        "N_AMBIGUOUS"          
# mdi-web-server-app-server-1              | [31] "N_DOWNSAMPLED"         "N_COLLAPSED"           "JXN_SEQ"              
# mdi-web-server-app-server-1              | [34] "MERGE_LEN"             "MICROHOM_LEN"          "JXN_BASES"            
# mdi-web-server-app-server-1              | [37] "FLANK_LEN1"            "FLANK_LEN2"            "N_CLUSTERED_JUNCTIONS"
# mdi-web-server-app-server-1              | [40] "SV_SIZE"               "GEN_REF_1"             "GEN_REF_2"            
# mdi-web-server-app-server-1              | [43] "GEN_COV_1"             "GEN_COV_2"             "TARGET_REGION"        
# mdi-web-server-app-server-1              | [46] "TARGET_POS_1"          "TARGET_POS_2"          "edgeType"             
# mdi-web-server-app-server-1              | [49] "chromIndex1"           "chromIndex2"           "strand1"              
# mdi-web-server-app-server-1              | [52] "strand2"               "insertSize"            "nSamples"             
# mdi-web-server-app-server-1              | [55] "nInstances"            "nSequenced"            "mapQ"                 
# mdi-web-server-app-server-1              | [58] "flankLength"           "nLinkedJunctions"      "samples"              
# mdi-web-server-app-server-1              | [61] "node1"                 "node2"                 "isCanonical"          
# mdi-web-server-app-server-1              | [64] "cChromIndex1"          "cChromIndex2"          "cRefPos1"             
# mdi-web-server-app-server-1              | [67] "cRefPos2"              "cStrand1"              "cStrand2"             
# mdi-web-server-app-server-1              | [70] "nodeCenter"            "refPosCenter"          "size"                 
# mdi-web-server-app-server-1              | [73] "cChrom1"               "cChrom2"               "X"                    
# mdi-web-server-app-server-1              | [76] "sourceId"              "sourceGenome"          "Project"              
# mdi-web-server-app-server-1              | [79] "Sample_ID"             "Clonal"                "SV_Type"              
# mdi-web-server-app-server-1              | [82] "genome1"               "genome2"               "Genome"               
# mdi-web-server-app-server-1              | [85] "group"

#----------------------------------------------------------------------
# junction expansion
#----------------------------------------------------------------------
expansionTableData <- reactiveVal(NULL)
expansionUIContents <- reactiveVal()
observeEvent(matchingJunctionsSelected(), {
    jxn <- matchingJunctionsSelected()
    svWGS_expandJunction_(
        jxn = jxn,
        targetId = jxn$sourceId,
        objectTableFn    = NULL, 
        expansionTableFn = expansionTableData, 
        expansionUIFn    = expansionUIContents,
        Alignment_Mode = "Evidence Consensus" # All Molecules  Reference Molecule
    )
    stopSpinner(session)
})
expansionTable <- bufferedTableServer(
    "expansionTable",
    id,
    input,
    tableData = expansionTableData,
    selection = 'none',
    options = list(
    )
)
output$expansionUI <- renderUI({
    expansionUIContents()
})
output$expansionSequence = renderUI({
    jxn <- matchingJunctionsSelected()
    tags$div(
        style = "border-top: 1px solid grey; margin-top: 10px; padding-top: 10px; width: 750px; word-wrap:break-word;",
        jxn$JXN_SEQ
    )
})
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    insertSizeDensityPlot$settings$replace(bm$outcomes$insertSizeDensityPlotSettings)
    sizeCorrelationPlot$settings$replace(bm$outcomes$sizeCorrelationPlotSettings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = list(
        insertSizeDensityPlotSettings = insertSizeDensityPlot$settings$all_,
        sizeCorrelationPlotSettings = sizeCorrelationPlot$settings$all_
    ),
    # isReady = reactive({ getStepReadiness(options$source, gcBiasModels()) }),
    # getGcBiasModels = getGcBiasModels_externalCall,
    # getBinNormalizedCN = getBinNormalizedCN,
    # getCnvJxnsNormalizedCN = getCnvJxnsNormalizedCN,
    # getNormalizedHmmCnvs = getNormalizedHmmCnvs,
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
