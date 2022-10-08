#----------------------------------------------------------------------
# server components for the sampleSummary appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
sampleSummaryServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'sampleSummary'
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

#----------------------------------------------------------------------
# sample selection and loading
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer(
    "source", 
    selection = "single"
)
projectName <- projectNameReactive(sourceId)
sampleData <- sampleDataReactive(sourceId)

#----------------------------------------------------------------------
# metric boxes (a series of info boxes)
#----------------------------------------------------------------------
metricBoxCss <- "
.metric-div {
    white-space: nowrap;
    display: inline-block;
    margin-right: 10px;
}
.metric-sub-div {
    height: 50px;
    display: inline-block;
    margin: 0 15px;
    float: left;    
}
.metric-icon {
    line-height: 50px;
    font-size: 1.5em;
}
.metric-contents p {
    line-height: 25px;
    margin: 0;
}
"
metricBox <- function(title, value, iconName, color) tags$div(
    class = "metric-div",
    tags$div(
        class = "metric-sub-div metric-icon",
        style = paste("color:", color),
        icon(iconName, verify_fa = FALSE)
    ),
    tags$div(
        class = "metric-div metric-contents",
        tags$p(tags$strong(title)),
        tags$p(value)        
    )
)
output$metrics <- renderUI({
    sampleData <- sampleData()
    req(sampleData)
    d <- sampleData$constants
    cellCount <- function(x) paste0(d[[x]], " (", as.integer(d[[x]] / d$num_cells * 100), "%)")
    tagList(
        tags$style(metricBoxCss),    
        fluidRow(
            style = "margin-bottom: 15px;",
            metricBox("Total Cells", d$num_cells, "circle", "blue"),
            metricBox("Good Cells",  cellCount("num_good_cells"), "check",  "green"),
            metricBox("Bad Cells",   cellCount("num_bad_cells"),  "close",  "#cc0000")
        )
    )
})

#----------------------------------------------------------------------
# aggregate plots of all good cells
#----------------------------------------------------------------------
windowSizes <- staticPlotBoxServer(
    "windowSizes",
    #----------------------------
    maxHeight = "400px",
    immediate = TRUE,
    #----------------------------
    envir = parent.frame(),
    settings = NULL,
    create = function(){
        sampleData <- sampleData()
        req(sampleData)
        d <- sampleData$colData[rejected == FALSE, .(N = .N), by = window_size]
        windowSizes$initializeFrame(
            title = projectName(),
            xlim = range(d$window_size, na.rm = TRUE),
            ylim = c(0, max(d$N, na.rm = TRUE) * 1.1),
            xlab = "Window Size (# of 20kb bins)",
            ylab = "# of Cells"
        )
        windowSizes$addPoints(
            x = d$window_size,
            y = d$N,
            typ = "h"
        )
    }
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

# red
# yellow
# aqua
# blue
# light-blue
# green
# navy
# teal
# olive
# lime
# orange
# fuchsia
# purple
# maroon
# black
