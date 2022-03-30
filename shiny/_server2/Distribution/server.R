
#----------------------------------------------------------------------
# server.R returns a function that defines reactive server actions
#----------------------------------------------------------------------

library(data.table)

# handle user actions and update output
server <- function(input, output, session){
    
    # set the user/session-specific environment
    sessionEnv <- environment()
    verbose <<- TRUE

    # declare session-specific data resources; only visible within server block
    # see: https://shiny.rstudio.com/articles/scoping.html
    initializeSession <- function(){
        reportProgress('initializeSession')
        sourceScripts(c(
            "../data_sources_server.R",
            "../data_sources_session.R",
            "../presets.R",
            "data_sources.R",
            "plot_data.R"
        ), env=sessionEnv)        
    }
    initializeSession()

    # load all samples and supporting data again when requested (e.g. when new data is available)
    observeEvent(input$refreshServer, {
        initalizeApp()
        initializeSession()
        updateSelectInput(session, 'project', choices=projects, selected='-')
    })    

    # clear the current filter selections
    observeEvent(input$resetFilters, {
        resetCommonFilters()
        #updateRadioButtons(session, 'plotColors', selected=names(plotColors)[1])
        #updateRadioButtons(session, 'pointCex', selected='1')
        #updateRadioButtons(session, 'jxnDisplayType', selected='Genome')
        #updateRadioButtons(session, 'charPerLine', selected=100)
    })
    
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
        
        print(loadCaptureTargets)
        
        loadCaptureTargets()
        #loadExpectedCnvs()
        #updateSelectInput(session, 'captureTarget', choices=captureTargetNames, selected='-')
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
        makeDistPlot()
    })
    output$eventCount <- renderUI({
        paste(eventCount_$SVs, 'SVs,',
              eventCount_$gapSplit, 'gaps+splits,',
              eventCount_$molecules, 'molecules,',
              eventCount_$alignments, 'alignments')
    })

}

