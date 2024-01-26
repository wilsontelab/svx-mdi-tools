#----------------------------------------------------------------------
# server components for the sampleExplorer appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svWGS_sampleExplorerServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'sampleExplorer'
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
    sourceSelection = "single",
    sampleGenomesSelection = "multiple"
)

#----------------------------------------------------------------------
# correlation plot of junction counts between samples
#----------------------------------------------------------------------
countCorrelationPlot <- staticPlotBoxServer(
    "countCorrelationPlot",
    margins = TRUE,
    title = TRUE,
    points = TRUE,
    legend = TRUE,
    settings = list(
        Limits = list(
            Min_Count = list(
                type = "numericInput",
                value = 0,
                min = 0, 
                max = 100,
                step = 10
            ),
            Max_Count = list(
                type = "numericInput",
                value = 100,
                min = 0, 
                max = 100,
                step = 10
            )
        )
    ), 
    size = "m",
    create = function() {
        sampleGenomes <- data$sampleGenomes()
        req(sampleGenomes)
        sampleIds <- sampleGenomes[, unique(Sample_ID)]
        req(length(sampleIds) == 2)
        jxns <- data$sampleGenomeJunctions()
        req(jxns, nrow(jxns) > 0)
        lim <- c(
            countCorrelationPlot$settings$get("Limits","Min_Count"),
            countCorrelationPlot$settings$get("Limits","Max_Count")
        )
        title <- countCorrelationPlot$settings$get("Plot_Frame", "Title", NULL)
        totalJxns <- paste(nrow(jxns), "Junctions")
        title <- if(is.null(title) || title == "") totalJxns else paste0(title, " (", totalJxns, ")")
        countCorrelationPlot$initializeFrame(
            xlim = lim,
            ylim = lim,
            xlab = sampleIds[1],
            ylab = sampleIds[2],
            # xaxs = "i",
            # yaxs = "i",
            title = title,
            cex.main = 0.95
        )
        countCorrelationPlot$addPoints(
            x = jitter(jxns[[sampleIds[1]]], a = 0.5),
            y = jitter(jxns[[sampleIds[2]]], a = 0.5),
            col = rgb(0, 0, 0, alpha = 0.2)
        )
        abline(0, 1)
        # countCorrelationPlot$addMarginLegend(
        #     xlim[2] * 1.1, maxY, lty = 1, lwd = 1.5, 
        #     legend = paste0(d$grouping$groups, " (", trimws(commify(d$grouping$groupCounts)), ")"),
        #     col = unlist(CONSTANTS$plotlyColors[1:d$grouping$nGroups]),
        #     bty = "n"
        # )
        stopSpinner(session)
    }
)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    countCorrelationPlot$settings$replace(bm$outcomes$countCorrelationSettings)
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
        countCorrelationSettings = countCorrelationPlot$settings$all_
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
