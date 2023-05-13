getSummaryAxisLabel <- function(column){
    if(column %in% c("eventSize","insertSize")){
        paste("log10(", column, ")")
    } else {
        column
    }
}
getSummaryAxisLines <- function(input){
    getLines <- function(column){
        if(column %in% c("eventSize","insertSize")){
            -10:10
        } else if(column %in% c("nInstances")) {
            seq(0, 100, 10)
        } else {
            0:1
        }
    }
    list(
        v = getLines(input$xAxisColumn),
        h = getLines(input$yAxisColumn),
        col = rgb(0, 0, 0, 0.15)
    )
}
summaryPlotServer <- function(plotNodes, input) mdiInteractivePlotServer(
    "summaryPlot",      
    click = TRUE,
    contents = reactive({ 
        d <- plotNodes()
        list(
            plotArgs = list(
                x = d$xj,
                y = d$yj,
                pch = 19,
                cex = 0.25,
                col = d$color,
                xlab = getSummaryAxisLabel(input$xAxisColumn),
                ylab = getSummaryAxisLabel(input$yAxisColumn),
                xaxs = "i",
                yaxs = "i"
            ),
            layout = list(
                width = 800,
                height = 800,
                pointsize = 8, # defaults to 8
                dpi = 96, # defaults to 96
                mai = c(0.75, 0.75, 0.1, 0.1),
                xlim = range(d$xj, na.rm = TRUE), # OR can be read from plotArgs
                ylim = range(d$yj, na.rm = TRUE)
            ),
            abline = getSummaryAxisLines(input) # optional named arguments passed to abline() after calling plot(plotArgs)
        ) 
    })
)

moleculePlotServer <- function(sourceId, segments, input, session) mdiInteractivePlotServer(
    "moleculePlot",      
    click = FALSE,
    #     parseLayout = function(x, y) list(x, y, layout) # to convert to plot space in a multi-plot layout
    contents = reactive({
        d <- segments()
        startSpinner(session, message = "building molecule plot")        
        pngFile <- file.path(sessionDirectory, "moleculePlotServer.png")
        png(
            pngFile,
            width = 800,
            height = 600,
            units = "px",
            pointsize = 8,
            res = 96,
            type = "cairo"
        )
        plotSegments(sourceId, d)
        dev.off()
        stopSpinner(session)
        list(
            pngFile = pngFile,
            layout = list(
                width = 800,
                height = 600,
                pointsize = 8, # defaults to 8
                dpi = 96, # defaults to 96
                mai = c(0.75, 0.75, 0.1, 0.1)
                # ,
                # xlim = range(d$xj, na.rm = TRUE), # OR can be read from plotArgs
                # ylim = range(d$yj, na.rm = TRUE)
            )
        ) 
    })
)
