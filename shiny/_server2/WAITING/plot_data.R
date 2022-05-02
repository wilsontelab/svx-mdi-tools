
#----------------------------------------------------------------------
# plot_data.R controls the data visualization actions
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# common plot values and colors
#----------------------------------------------------------------------
getColors <- function(){  
    switch(input$plotColors,
        duplex = {
            ifelse(sv$N_DUPLEX > 0, plotColors$duplex$duplex, plotColors$duplex$single)
        },
        svType = sapply(sv$JXN_TYPE, function(jt){
            jt_ <- jxnTypes[[jt]]
            if(is.null(jt_)) return(plotColors$svType$Undet)
            c <- plotColors$svType[[jt_]]
            if(is.null(c)) return(plotColors$svType$Undet)
            c
        }),
        sample = sapply(sv$SAMPLE, function(smp){
            for(i in 1:nSamples) {
                if(smp == input[[paste0('sample_', i)]]) return(plotColors$sample[i])
            }
            return("grey")
        })
    ) 
}

#----------------------------------------------------------------------
# function to wrap a DT table to add some desired formatting
#----------------------------------------------------------------------
formatDTTable <- function(table, toCommify, pageLength=5, lengthMenu=c(5,10,25,50,100)){
    reportProgress('formatDTTable')
    if(!is.null(table)) formatRound(
        datatable(
            table,
            options=list(
                pageLength=pageLength,
                lengthMenu=lengthMenu,
                searching=TRUE
            ),
            selection = 'single',
            rownames= FALSE                       
        ),
        toCommify,
        digits = 0
    ) else table    
}
#----------------------------------------------------------------------
# functions for each of the major rows of the data output
#----------------------------------------------------------------------

# multiple plot views of groups of SVs (top row of display)
svPlot <- NA # the SV data set from which plot will select subsets
sv     <- NA # the working SV data set for each plot sub
makeSVPlots <- function(){
    reportProgress('makeSVPlots')
    # get the requested SVs to plot
    if(input$sample_1 == '-' | input$sample_1 == '') return()
    svPlot <<- getFilteredSVs(c(
        'center','size', # for triangle plot
        'j1','j2','ct1','ct2', # for linear plot
        'MICROHOM_LEN', 'FRAC_SHARED_PROPER', # for microhomology plot
        'SAMPLE' # for point coloring
    ))
    
    # create a horizontal layout of all requested plot types
    nPlots <- length(input$showPlots)
    if(nPlots == 0) return()
    layout(matrix(1:nPlots, nrow=1, byrow=TRUE))
    par(mar=c(5.1,5.1,0.1,0.1))
    if('triangle' %in% input$showPlots) makeTrianglePlot()
    if('linear' %in% input$showPlots) makeLinearPlot()    
    if('corr' %in% input$showPlots) makeCorrelationPlot()
    if('size_uHom' %in% input$showPlots) makeSizeUHomPlot()
    if('uHom' %in% input$showPlots) makeMicrohomologyHistogram()    
}

# one table of filtered SVs 
svTable <- NA # SV data set from which table will be assembled and juncions will be recovered
makeSVTable <- function(){
    reportProgress('makeSVTable')
    if(input$sample_1 == '-' | input$sample_1 == '') return()
    svTable <<- getFilteredSVs(names(svTableNames))
    #svTable$SV_SIZE <<- ifelse( # calculate SV_SIZE for gaps (not currently done by find)
    #    svTable$SV_SIZE != 0 | svTable$SV_SIZE == 'T',
    #    svTable$SV_SIZE,
    #    abs(svTable$PROX_JXN_POS2 - svTable$PROX_JXN_POS1)
    #)
    svTable <<- svTable[,names(svTableNames),with=FALSE] # remove columns used to filter but not displayed
    colnames(svTable) <<- unname(unlist(svTableNames)) # override column names in table for display
    svTable # let Shiny DT build the final table
}



# plot of fragment sizes and junction positions to judge evidence uniqueness 
makeZoomPlot <- function(){
    reportProgress('makeZoomPlot')
    zoom <- zoomInfo()
    if(is.null(zoom)) return()
    #par(mar=c(4.1,4.1,1.1,1.1))
    #plotMolSizePos(zoom)
    NULL
}

# alignment view of a single SV junction
nullImage <- list(src="")
jxnScalars <- 1:10
maxJxnWidth <- 1600
makeJunctionConsensus <- function(input){
    reportProgress('makeJunctionConsensus')
    jxnImgData <- jxnImgData()
    zoom <- zoomInfo() 
    if(is.null(jxnImgData) || is.null(zoom)) return(nullImage)    
    png <- getSVJunctionImage(input$project, zoom$smp, zoom$svId, 'consensus')
    dim <- createJunctionConsensus(jxnImgData, zoom, png, input$clipMode)
    xScalar <- max(1, which(jxnScalars * dim[1] <= maxJxnWidth))    
    list(src=png, width=xScalar * dim[1], height=dim[2] * 6, style="image-rendering: pixelated;")
}
makeJunctionZoom <- function(input){
    reportProgress('makeJunctionZoom')
    zoom <- zoomInfo() 
    if(is.null(zoom) || !isSequencedJunction(zoom$sv)) {
        jxnImgData(NULL)
        return(nullImage)
    }
    png <- getSVJunctionImage(input$project, zoom$smp, zoom$svId)
    dim <- createJunctionZoom(zoom, png, input$clipMode)
    xScalar <- max(1, which(jxnScalars * dim[1] <= maxJxnWidth))    
    yScalar <- max(1, which(jxnScalars * dim[2] <= 200))
    list(src=png, width=xScalar * dim[1], height=yScalar * dim[2], style="image-rendering: pixelated;")
}

# one table of molecules for a selected SV
makeNodeTable <- function(){
    reportProgress('makeNodeTable')
    if(input$sample_1 == '-') return()
    zoom <- zoomInfo()
    if(is.null(zoom)) return()
    inNames <- names(nodeTableNames)
    zoom$nodes <- zoom$nodes[IS_JUNCTION_NODE==TRUE,..inNames]
    setnames(zoom$nodes, unname(unlist(nodeTableNames)))  
    zoom$nodes # let Shiny DT build the final table
}

