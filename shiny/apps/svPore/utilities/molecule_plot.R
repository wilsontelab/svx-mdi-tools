#----------------------------------------------------------------------
# dot plot of a single molecule's alignments
#----------------------------------------------------------------------
# plotting constants
dotPlotColors <- list(
    M = CONSTANTS$plotlyColors$grey,
    m = CONSTANTS$plotlyColors$black,
    I = CONSTANTS$plotlyColors$red,
    D = CONSTANTS$plotlyColors$blue,
    H = CONSTANTS$plotlyColors$orange,
    S = NA
)
dotStackOrder <- list(
    M = 1,
    m = 2,
    I = 3,
    D = 4,
    H = 5,
    S = 6
)

# use CIGAR string and refPos to parse the x-y coordinates of an alignment dot plot
cigarDotPlot <- function(cigar, qryPos, refPos, strand){
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    nOps <- length(operations)
    dt <- do.call(rbind, lapply(1:nOps, function(i){
        length <- lengths[i]
        operation <- operations[i]

        # deleted segments are present in reference, not in query
        if(operation == 'D'){
            dt <- data.table(
                x = qryPos,                
                y = refPos:(refPos + length - 1),
                operation = operation,
                level = 1
            )
            refPos <<- refPos + length
            
        # inserted bases are present in query, but not reference
        } else if(operation == 'I' || operation == 'S'){
            dt <- data.table(
                x = qryPos:(qryPos + length - 1),                
                y = refPos,
                operation = if(operation == 'S') operation else operation,
                level = 3
            )
            qryPos <<- qryPos + (if(operation == 'S') 0 else length)
            
        # M or m operations present in reference and query
        } else {
            dt <- data.table(
                x = qryPos:(qryPos + length - 1),                
                y = refPos:(refPos + length - 1),
                operation = "M",
                level = 2
            )
            refPos <<- refPos + length
            qryPos <<- qryPos + length
        }
        dt
    }))
    notS <- dt$operation != "S"
    if(strand < 0){
        nNotS <- sum(notS)
        dt$x[notS] <- dt$x[notS][nNotS:1]
    }
    list(
        qryPos = qryPos,
        refPos = refPos,
        dt = dt,
        x = mean(dt$x[notS]),
        y = mean(dt$y[notS])
    )
}

# build one locus plot points without plotting yet, may reorder loci later
buildLocus <- function(segments, minQStarts){
    pts <- data.table(x = integer(), y = integer(), col = character(), 
                      operation = character(), type = character())
    lbls <- data.table(x = integer(), y = integer(), label = character())
    ylim <- range(segments[, rStart], segments[, rEnd], na.rm = TRUE)
    bwls <- data.table(intercept = integer(), slope = integer())
    tryCatch({
        for(qName_ in segments[, unique(qName)]){
            xq <- segments[qName == qName_]
            minQs <- minQStarts[qName == qName_, qStart]
            bwls <- rbind(bwls, xq[1, .(
                intercept = if(strand > 0) min(rStart, rEnd) - qStart + minQs 
                                      else max(rStart, rEnd) + qStart - minQs, 
                slope = strand
            )])

            for(edge_ in xq[, edge]){
                xe <- xq[edge == edge_]
                if(xe$edgeType == edgeTypes$ALIGNMENT){
                    qryPos <- xe[, qStart]
                    refPos <- xe[, min(rStart, rEnd)]
                    dp <- cigarDotPlot(xe[, cigar], qryPos, refPos, xe[, strand])
                    pts <- rbind(pts, data.table(
                        x = dp$dt$x,
                        y = dp$dt$y,
                        col = unlist(dotPlotColors[dp$dt$operation]),
                        operation = dp$dt$operation,
                        type = "A"
                    ))
                    lbls <- rbind(lbls, data.table(
                        x = mean(dp$dt$x), 
                        y = ylim[2], 
                        label = xe[, label]
                    ))

                } else {
                    qryPos <- xe[, qStart]
                    refPos <- xe[, rStart]
                    endQryPos <- qryPos + xe[, qOffset - 1]
                    operation <- xe[, if(qOffset > 1) "I" else "H"]
                    pts <- rbind(pts, data.table(
                        x = qryPos:endQryPos,
                        y = refPos,
                        col = unlist(dotPlotColors[operation]),
                        operation = operation,
                        type = "J"
                    ))
                    pts <- rbind(pts, data.table(
                        x = endQryPos,
                        y = xe[, rStart]:xe[, rEnd],
                        col = dotPlotColors$D,
                        operation = "D",
                        type = "J"
                    ))                  
                }
            }        
        }
    }, error = function(e){
        # str(xe)
        print(e)
    })
    list(
        chrom = segments[1, chrom],
        ylim = ylim,  
        points = pts[order(unlist(dotStackOrder[pts$operation]))],
        labels = lbls,
        bandwidthLines = bwls    
    )   
}

# add interlocus segments to the locus plots
buildInterLocus <- function(segments, plotData){
    N <- nrow(segments)
    if(N == 0) return(plotData)
    tryCatch({
        for(i in 1:N){
            x <- segments[i]
            isUp <- x[, if(locus2 > locus1) TRUE else FALSE]

            locus <- as.character(x[, locus1])
            pts <- plotData[[locus]]$points
            qryPos <- x[, qStart]
            endQryPos <- qryPos + x[, qOffset - 1]
            refPos <- x[, rStart]
            limitPos <- plotData[[locus]]$ylim[if(isUp) 2 else 1]
            operation <- x[, if(qOffset > 1) "I" else "H"]
            plotData[[locus]]$points <- rbind(plotData[[locus]]$points, data.table(
                x = qryPos:endQryPos,
                y = refPos,
                col = unlist(dotPlotColors[operation]),
                operation = operation,
                type = "J"
            ))
            plotData[[locus]]$points <- rbind(plotData[[locus]]$points, data.table(
                x = endQryPos,
                y = refPos:limitPos,
                col = dotPlotColors$D,
                operation = "D",
                type = "J"
            )) 

            locus <- as.character(x[, locus2])
            pts <- plotData[[locus]]$points
            refPos <- x[, rEnd]
            limitPos <- plotData[[locus]]$ylim[if(isUp) 1 else 2]
            plotData[[locus]]$points <- rbind(plotData[[locus]]$points, data.table(
                x = endQryPos,
                y = refPos:limitPos,
                col = dotPlotColors$D,
                operation = "D",
                type = "J"
            ))                 
        }
    }, error = function(e){
        str(segments)
        str(plotData)
        print(e)
    })
    plotData
}

# plot a single chromosome segment
renderAlignmentPlot <- function(xlim, d){  
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = xlim,
        ylim = d$ylim,
        xlab = "",
        ylab = d$chrom,
        xaxt = "n"
    )  
    for(i in 1:nrow(d$bandwidthLines)) abline(
        d$bandwidthLines[i, intercept], 
        d$bandwidthLines[i, slope],
        col = CONSTANTS$plotlyColors$green,
        lty = 2
    )    
    points(d$points$x, d$points$y, col = d$points$col, pch = 19, cex = 0.25)
    for(i in 1:nrow(d$labels)) d$labels[i, text(x, y, label)]
}

# plot the read(s) QUAL profile(s)
renderBaseQualPlot <- function(xlim, read){
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = xlim,
        ylim = c(0, 45),
        xlab = "Position in Read",
        ylab = "QUAL",
        xaxt = "n"
    )
    abline(h = seq(10, 40, 10), col = CONSTANTS$plotlyColors$grey)
    x <- xlim[1]:xlim[2]
    points(
        x = x,
        y = read$QUAL[x],
        pch = 19,
        cex = 0.25,
        col = rgb(0,0,0,0.1)
    )
}

# plot the read(s) M/D/I profile(s)
renderReadQualPlot <- function(xlim, segments){
    par(mar = c(4.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = xlim,
        ylim = c(0, 4),
        xlab = "Position in Read",
        ylab = "Alignment"
    )
    abline(h = 1:3, col = CONSTANTS$plotlyColors$grey)
    xq <- segments$dt[qName == segments$qNames[1]]
    dt <- do.call(rbind, lapply(xq[edgeType == edgeTypes$ALIGNMENT, edge], function(edge_){
        xe <- xq[edge == edge_]
        qryPos <- xe[, qStart]
        refPos <- xe[, min(rStart, rEnd)]
        dp <- cigarDotPlot(xe[, cigar], qryPos, refPos, xe[, strand])
        dp$dt[, .(x = x, y = level)]        
    }))
    colors <- c(CONSTANTS$plotlyColors$blue, CONSTANTS$plotlyColors$black, CONSTANTS$plotlyColors$red)
    points(
        x = dt$x,
        y = jitter(dt$y, amount = 0.4),
        pch = 19,
        cex = 0.25,
        col = colors[dt$y]
    )
}

# render the composite QC plot; this is the main function call, cascading upwards
baseQualHeight <- 0.1
readQualHeight <- 0.19
plotSegments <- function(sourceId, segments) {
    xlim <- range(segments$dt[, qStart], segments$dt[, qEnd], na.rm = TRUE)

    # dprint(segments$qNames)

    # prepare for alignment to multiple genome windows, i.e., loci

    qualRowHeightSum <- baseQualHeight + readQualHeight
    alnRowHeight <- (1 - qualRowHeightSum) / segments$nLoci
    heights <- c(rep(alnRowHeight, segments$nLoci), baseQualHeight, readQualHeight)
    nPlotRows <- segments$nLoci + 2    
    layout(matrix(1:nPlotRows, ncol = 1), heights = heights) 

    # plot alignments
    minQStarts <- segments$dt[edgeType == edgeTypes$ALIGNMENT, .(qStart = min(qStart, na.rm = TRUE)), by = .(qName)]
    plotData <- lapply(segments$uniqueLoci, function(locus_) buildLocus(segments$dt[locus1 == locus_ & locus2 == locus_], minQStarts) )
    names(plotData) <- as.character(segments$uniqueLoci)
    orderedLoci <- rev(sort(segments$uniqueLoci))

    if(segments$nLoci > 1) {
        pairs <- as.data.table(expand.grid(leftLocus = segments$uniqueLoci, rightLocus = segments$uniqueLoci))
        pairs <- pairs[leftLocus != rightLocus]
        apply(pairs, 1, function(pair){
            plotData <<- buildInterLocus(segments$dt[locus1 == pair[1] & locus2 == pair[2]], plotData)
        })
    }
    # for(i in 1:segments$maxI) if(segments$dt[i, locus1 != locus2]){
    #     leftLocus  <- segments$dt[i, locus1]
    #     rightLocus <- segments$dt[i, locus2] 
    # }
    for(locus_ in orderedLoci) renderAlignmentPlot(xlim, plotData[[as.character(locus_)]])

    # load read data from indexed file    
    read <- getReadFromSequenceIndex(sourceId, segments$qNames[1])

    # plot QUAL
    renderBaseQualPlot(xlim, read)
    renderReadQualPlot(xlim, segments)
}
