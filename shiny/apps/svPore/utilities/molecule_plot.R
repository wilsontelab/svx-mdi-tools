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
buildLocus <- function(md, minQStarts){
    pts <- data.table(x = integer(), y = integer(), col = character(), 
                      operation = character(), type = character())
    lbls <- data.table(x = integer(), y = integer(), label = character())
    ylim <- range(md[, rStart], md[, rEnd], na.rm = TRUE)
    bwls <- data.table(intercept = integer(), slope = integer())
    tryCatch({
        for(readKey_ in md[, unique(readKey)]){
            xq <- md[readKey == readKey_]
            minQs <- minQStarts[readKey == readKey_, qStart]
            bwls <- rbind(bwls, xq[1, .(
                intercept = if(strand > 0) min(rStart, rEnd) - qStart + minQs 
                                      else max(rStart, rEnd) + qStart - minQs, 
                slope = strand
            )])

            for(edgeN_ in xq[, edgeN]){
                xe <- xq[edgeN == edgeN_]
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
        chrom = md[1, chrom],
        ylim = ylim,  
        points = pts[order(unlist(dotStackOrder[pts$operation]))],
        labels = lbls,
        bandwidthLines = bwls    
    )   
}

# add interlocus segments to the locus plots
buildInterLocus <- function(md, pd){
    N <- nrow(md)
    if(N == 0) return(pd)
    tryCatch({
        for(i in 1:N){
            x <- md[i]
            isUp <- x[, if(locus2 > locus1) TRUE else FALSE]

            locus <- as.character(x[, locus1])
            pts <- pd[[locus]]$points
            qryPos <- x[, qStart]
            endQryPos <- qryPos + x[, qOffset - 1]
            refPos <- x[, rStart]
            limitPos <- pd[[locus]]$ylim[if(isUp) 2 else 1]
            operation <- x[, if(qOffset > 1) "I" else "H"]
            pd[[locus]]$points <- rbind(pd[[locus]]$points, data.table(
                x = qryPos:endQryPos,
                y = refPos,
                col = unlist(dotPlotColors[operation]),
                operation = operation,
                type = "J"
            ))
            pd[[locus]]$points <- rbind(pd[[locus]]$points, data.table(
                x = endQryPos,
                y = refPos:limitPos,
                col = dotPlotColors$D,
                operation = "D",
                type = "J"
            )) 

            locus <- as.character(x[, locus2])
            pts <- pd[[locus]]$points
            refPos <- x[, rEnd]
            limitPos <- pd[[locus]]$ylim[if(isUp) 1 else 2]
            pd[[locus]]$points <- rbind(pd[[locus]]$points, data.table(
                x = endQryPos,
                y = refPos:limitPos,
                col = dotPlotColors$D,
                operation = "D",
                type = "J"
            ))                 
        }
    }, error = function(e){
        str(md)
        str(pd)
        print(e)
    })
    pd
}

# plot a single chromosome segment
renderAlignmentPlot <- function(xlim, pd){ # pd = plotData for one locus
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = xlim,
        ylim = pd$ylim,
        xlab = "",
        ylab = pd$chrom,
        xaxt = "n",
        xaxs = "i", yaxs = "i"
    )  
    tryCatch({
        for(i in 1:nrow(pd$bandwidthLines)) abline(
            pd$bandwidthLines[i, intercept], 
            pd$bandwidthLines[i, slope],
            col = CONSTANTS$plotlyColors$green,
            lty = 2
        )          
    }, error = function(e) NULL)
    points(pd$points$x, pd$points$y, col = pd$points$col, pch = 19, cex = 0.25)
    # for(i in 1:nrow(pd$labels)) pd$labels[i, text(x, y, label)]
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
renderReadQualPlot <- function(xlim, md){
    par(mar = c(4.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = xlim,
        ylim = c(0, 4),
        xlab = "Position in Read",
        ylab = "Alignment"
    )
    abline(h = 1:3, col = CONSTANTS$plotlyColors$grey)
    xq <- md$dt[readKey == md$readKey[1]]
    dt <- do.call(rbind, lapply(xq[edgeType == edgeTypes$ALIGNMENT, edgeN], function(edgeN_){
        xe <- xq[edgeN == edgeN_]
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
renderMoleculePlot <- function(md) { # md == molecule(s) data
    xlim <- range(md$dt[, qStart], md$dt[, qEnd], na.rm = TRUE)

    # prepare for alignment to multiple genome windows, i.e., loci
    # qualRowHeightSum <- baseQualHeight + readQualHeight
    # alnRowHeight <- (1 - qualRowHeightSum) / md$nLoci
    # heights <- c(rep(alnRowHeight, md$nLoci), baseQualHeight, readQualHeight)
    # nPlotRows <- md$nLoci + 2    
    # layout(matrix(1:nPlotRows, ncol = 1), heights = heights) 

    # alnRowHeight <- (1 - readQualHeight) / md$nLoci
    # heights <- c(rep(alnRowHeight, md$nLoci), readQualHeight)
    # nPlotRows <- md$nLoci + 1   
    # layout(matrix(1:nPlotRows, ncol = 1), heights = heights) 

    alnRowHeight <- 1 / md$nLoci
    heights <- rep(alnRowHeight, md$nLoci)
    nPlotRows <- md$nLoci
    layout(matrix(1:nPlotRows, ncol = 1), heights = heights) 

    # plot alignments
    startSpinner(session, message = "calculating plot data")
    minQStarts <- md$dt[edgeType == edgeTypes$ALIGNMENT, .(qStart = min(qStart, na.rm = TRUE)), by = .(readKey)]
    plotData <- lapply64(md$uniqueLoci, function(locus_) buildLocus(md$dt[locus1 == locus_ & locus2 == locus_], minQStarts) )
    names(plotData) <- as.character(md$uniqueLoci)
    orderedLoci <- rev(sort(md$uniqueLoci))
    if(md$nLoci > 1) {
        pairs <- as.data.table(expand.grid(leftLocus = md$uniqueLoci, rightLocus = md$uniqueLoci))
        pairs <- pairs[leftLocus != rightLocus]
        mapply(function(pair){
            plotData <<- buildInterLocus(md$dt[locus1 == pair[1] & locus2 == pair[2]], plotData)
        }, pairs)
    }

    startSpinner(session, message = "rendering plots")
    invisible(sapply64(orderedLoci, function(locus_) renderAlignmentPlot(xlim, plotData[[as.character(locus_)]])))

    # load read data from indexed file    
    # read <- getReadFromSequenceIndex(sourceId, md$qNames[1])

    # plot QUAL
    # renderBaseQualPlot(xlim, read)
    # renderReadQualPlot(xlim, md)
}
