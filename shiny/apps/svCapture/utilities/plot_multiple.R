
#----------------------------------------------------------------------
# plot_multiple.R creates different kinds of multi-SV representation plots
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# triangle plot of all padded target regions together
#----------------------------------------------------------------------
makeTrianglePlot <- function(){
    reportProgress('makeTrianglePlot')
    
    # initialize plot (always create plot, even if contains no SVs)
    axisLabelCex   <- 0.85
    regionLabelCex <- 0.85
    hasInterRegion <- FALSE
    for(tgtType in c('tt','ta','t-')) hasInterRegion <- hasInterRegion | tgtType %in% input$targetTypeFilter
    ymax <- if(hasInterRegion) maxI else 3 * projectInfo$REGION_PADDING
    plot(0, 0, type="n", bty="n", 
         xlab="", ylab="", xaxt="n", yaxt="n",
         xlim=c(1,maxI), ylim=c(1,ymax),
         xaxs="i", yaxs="i")    

    # add region labels, region-separating gridlines, and x-coordinate axes
    addGridline <- function(i, slope, col="black"){
        if(slope == 1){
            lines(c(i, i + (maxI - i) / 2), c(0, maxI - i), col=col)
        } else {
            lines(c(i, i / 2), c(0, maxI - (maxI - i)), col=col)
        }
    }
    lines(c(0,maxI),c(1,1), col="black", lwd=2)
    for(i in 1:nCaptureTargets){
        # diagonal lines
        i0 <- captureTargets[i,'iOffset']
        i1 <- i0 + projectInfo$REGION_PADDING
        i2 <- i1 + captureTargets[i,'length']
        i3 <- i2 + projectInfo$REGION_PADDING
        addGridline(i1, 1, col="grey60")
        addGridline(i1, -1, col="grey60")
        addGridline(i2, 1, col="grey60")
        addGridline(i2, -1, col="grey60")  
        addGridline(i0, 1)
        addGridline(i3, -1)    
        
        # coordinate axis per region
        unit <- 50000
        at <- seq(unit, captureTargets[i,'paddedLength']-unit, unit)
        axis(1, at=at+captureTargets[i,'iOffset'], labels=round((captureTargets[i,'start']+at)/1e6,2))
        
        # region name label
        mtext(paste(captureTargets[i,'region'], captureTargets[i,'chrom'], sep=","),
              side = 1, line = 2.25, at = captureTargets[i,'iOffset'] + captureTargets[i,'paddedLength'] / 2,
              cex=regionLabelCex)
        
        # coordinate name
        mtext("SV center (Mbp)", side = 1, line = 3.5, at = maxI / 2,
              cex=axisLabelCex)
    }

    # add circle around where we expect to see CNVs in the experiment
    points(expectedCnvs$center, expectedCnvs$size, pch=1,
           col="blue", cex=2)

    # and highlight the one select SV position in red
    zoom <- zoomInfo()
    if(!is.null(zoom)){
        points(zoom$sv$center, zoom$sv$size, pch=1, col="red", cex=2)
    }

    # add size axes to the y axis (to judge size of event)
    unit <- 50000
    at <- seq(unit, max(captureTargets$paddedLength, na.rm=TRUE)-unit, unit)
    axis(2, at=at, labels=round(at/1e3,2))
    mtext('SV size (kbp)', side = 2, line = 2.5, at = max(captureTargets$paddedLength)/2,
          cex=axisLabelCex)

    # add a legend
    if(input$plotColors == 'sample') {
        legend <- sapply(1:nSamples, function(i) input[[paste0('sample_', i)]])
        legend <- legend[legend != "-"]
        col <- plotColors$sample[1:length(legend)]
    } else {
        legend <- names(plotColors[[input$plotColors]])
        col <- unlist(plotColors[[input$plotColors]])
    }
    legend("topright", legend, col=col, pch=19, cex=1)    

    # ensure that both SV ends are in a padded region (not necessarily the same region)
    if(nrow(svPlot)==0) return()
    sv <<- svPlot[!is.na(svPlot$size),] 
    if(nrow(sv)==0) return();

    # add one data point per SV junction
    colors <- getColors()
    pch <- 16 # 16 18 sv$JXN_TYPE
    points(sv$center, sv$size, pch=pch, cex=as.numeric(input$pointCex) * 0.95, col=colors, font=2)
}

#----------------------------------------------------------------------
# linear plot of all padded target regions together, stacked on same x-axis
#----------------------------------------------------------------------
makeLinearPlot <- function(){
    reportProgress('makeLinearPlot')
    
    # ensure that both SV ends are in the SAME padded region
    if(nrow(svPlot)==0) return()
    sv <<- svPlot[!is.na(svPlot$j1) & !is.na(svPlot$j2) & svPlot$ct1==svPlot$ct2,]
    if(nrow(sv)==0) return();
    
    # reorder the SV ends into low>>high coordinate pairs on chromosome
    sv$j1_ <- pmin(sv$j1,sv$j2)
    sv$j2_ <- pmax(sv$j1,sv$j2)    
    sv <- sv[order(sv$j1_),]

    # stack the SV bars
    levelEnds <- rep(-1e10,10000)
    spacer <- 3*projectInfo$REGION_PADDING / 100
    coord <- mapply(function(j1, j2){
        stackLevel <- 1
        while(j1 <= levelEnds[stackLevel] + spacer){
            stackLevel <- stackLevel + 1
        }
        levelEnds[stackLevel] <<- max(levelEnds[stackLevel], j2)        
        c(j1,j2,stackLevel,stackLevel)
    }, sv$j1_, sv$j2_)

    # make the stacked SV bar plot
    colors <- getColors()
    matplot(coord[1:2,]/1e3, coord[3:4,],
        typ="l", lty=1, lwd=3, col=colors,
        bty="n", 
        xlab="Relative Position (kbp)", ylab="",
        xlim=c(1,3*projectInfo$REGION_PADDING)/1e3, ylim=c(0,max(which(levelEnds>0)+1)),
        xaxs="i", yaxs="i"
    )
    
    # add lines to denote boundaries of capture and padded regions
    abline(v=projectInfo$REGION_PADDING*0:3/1e3)
}

#----------------------------------------------------------------------
# simple histogram of the microhomology lengths (neg = insertion) of all sequenced junctions
#----------------------------------------------------------------------
makeMicrohomologyHistogram <- function(){
    reportProgress('makeMicrohomologyHistogram')
    
    # ensure that we only plot sequenced junctions (where we know microhomology info)
    if(nrow(svPlot)==0) return()
    sv <<- svPlot[svPlot$JXN_BASES != '0',]
    if(nrow(sv)==0) return()
    if(sum(!is.na(sv$MICROHOM_LEN)) == 0) return() # not clear why there is one NA row...
    
    # make histogram (ensure that we always have a bar for every length
    len <- aggregate(sv$MICROHOM_LEN, by=list(sv$MICROHOM_LEN), length)
    plot(0, 0, type="n", 
         xlim=uHomLim, ylim=c(0,max(len[[2]])*1.1), xlab=uHomLab, ylab="# of junctions")
    abline(v=0,col="grey")    
    points(len[[1]], len[[2]], type="h", col="green4", lwd=1.5)
}

#----------------------------------------------------------------------
# microhomology lengths (neg = insertion) of all sequenced junctions against
# the fraction of all molecule endpoints that matched a corresponding proper fragment
#---------------------------------------------------------------------- 
makeCorrelationPlot <- function(){
    reportProgress('makeCorrelationPlot')
    
    # ensure that we only plot sequenced junctions (where we know microhomology info)
    if(nrow(svPlot)==0) return()
    sv <<- svPlot[svPlot$JXN_BASES != '0',]
    if(nrow(sv)==0) return()         
    
    # make the correlation plot, jitter to help points be visible
    colors <- getColors()
    plot(jitter(sv$MICROHOM_LEN), jitter(sv$FRAC_SHARED_PROPER,amount=0.1), typ="p",
         xlim=uHomLim, ylim=c(-0.1,2.1),xlab=uHomLab,ylab="Avg. Proper Ends",
         pch=16, col=colors, cex=as.numeric(input$pointCex) * 1)
    abline(v=0)
    abline(h=0:2)    
}
makeSizeUHomPlot <- function(){
    reportProgress('makeSizeUHomPlot')
    
    # ensure that we only plot sequenced junctions (where we know microhomology info)
    if(nrow(svPlot)==0) return()
    sv <<- svPlot[svPlot$JXN_BASES != '0',]
    if(nrow(sv)==0) return()
    
    # make the correlation plot, jitter to help points be visible
    colors <- getColors()
    plot(jitter(sv$MICROHOM_LEN), log10(sv$size), typ="p",
         xlim=uHomLim, ylim=c(2,6), xlab=uHomLab, ylab="log10 SV Size",
         pch=16, col=colors, cex=as.numeric(input$pointCex) * 1)
    abline(v=0)
    abline(h=2:6)    
}

