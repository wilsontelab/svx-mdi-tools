# module server function, in myModule_server.R
myModuleServer <- function(id, ...) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------
    
    # clear the current filter selections
    observeEvent(input$resetFilters, { resetCommonFilters() })
    
    # update the project selection
    clearSamples <- function(){
        lapply(1:nSamples, function(i){
            updateSelectInput(session,
                              paste0('sample_', i),
                              choices = samples[[input$project]],
                              selected = '-')
        })
    }    
    observeEvent(input$project, {
        loadProjectInfo(input$project)
        loadCaptureTargets()
        loadExpectedCnvs()
        updateSelectInput(session, 'captureTarget', choices=captureTargetNames, selected='-')
        clearSamples()
    })
    
    # clear all current sample selections
    observeEvent(input$clearSamples, { clearSamples() })
    
    ## enable the Excel download of the current filtered data set
    #output$downloadExcel <- downloadHandler(
    #    filename = function() {
    #        paste("SV_data_", Sys.Date(), ".xlsx", sep="")
    #    },
    #    content = function(file) {
    #        svt <- getFilteredSVs(svExcelColumns)
    #        rows <- order(svt$SV_SIZE, svt$RNAME1, svt$PROX_JXN_POS1, svt$RNAME2, svt$PROX_JXN_POS2)
    #        write_xlsx(svt[rows,svExcelColumns], path=file, col_names=TRUE, format_headers=TRUE)  
    #    },
    #    contentType = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    #)
    
    # activate the row of multi-SV output plots
    output$plotOutput <- renderPlot({
        makeSVPlots()
    })
    
    # activate the table of multiple SVs
    output$svTableOutput <- renderDataTable(formatDTTable(
        makeSVTable(), svToCommify
    ))
    
    output$consensusPlot <- renderImage({
        makeJunctionConsensus(input)
    })
    output$junctionPlot <- renderImage({
        makeJunctionZoom(input)
    })
    
    ## activate the output rows dependent on a single selected junction
    #output$zoomPlots <- renderPlot({
    #    makeZoomPlot()
    #})
    #output$junctionZoom <- renderText(
    #    makeJunctionZoom()
    #)
    output$nodeTableOutput <- renderDataTable(formatDTTable(
        makeNodeTable(), nodeToCommify, pageLength=25
    ))

# module return value (if none, use NULL instead of a list)
NULL
# list(
#     myValue = reactiveVal(input$inputName),
#     myMethod = function(...) {}
# )
#----------------------------------------------------------------------
})}
