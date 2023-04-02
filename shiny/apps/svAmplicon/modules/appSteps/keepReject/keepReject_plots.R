#----------------------------------------------------------------------
# construct a composite paired quality plot, i.e. MAPQ vs. QUAL
#----------------------------------------------------------------------
pairedQualityPlotServer <- function(pathClassMoleculeTypes) mdiInteractivePlotServer(
    "pairedQualityPlot",       
    click = TRUE,
    contents = reactive({ 
        mts <- pathClassMoleculeTypes()        
        thresholds <- app$keepReject$thresholds()[[1]]
        list(
            plotArgs = list(
                jitter(mts$minBaseQual, amount = 1),
                jitter(mts$minMapQ,     amount = 1),
                pch = 19,  
                cex = 0.5,
                xlab = "min(QUAL)",
                ylab = "min(MAPQ)",
                xlim = c(0, 45),
                ylim = c(0, 65),
                col = ifelse(mts$passedQualityFilters, CONSTANTS$plotlyColors$green, CONSTANTS$plotlyColors$red)
            ),
            abline = list(h = thresholds$minMapQ, v = thresholds$minBaseQual),
            layout = list(
                width = 395,
                height = 395,
                pointsize = 9, # defaults to 8
                mai = c(0.75, 0.75, 0.1, 0.1)
            )
        ) 
    })
)

#----------------------------------------------------------------------
# quality profile and dot plot of a single molecule's reads and alignments
#----------------------------------------------------------------------

# plotting constants
dotPlotColors <- list(
    M = CONSTANTS$plotlyColors$grey,
    m = CONSTANTS$plotlyColors$orange,
    I = CONSTANTS$plotlyColors$red,
    D = CONSTANTS$plotlyColors$blue,
    H = CONSTANTS$plotlyColors$green,
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
                operation = "M" # operation
            )
            refPos <<- refPos + length
            
        # inserted bases are present in query, but not reference
        } else if(operation == 'I' || operation == 'S'){
            dt <- data.table(
                x = qryPos:(qryPos + length - 1),                
                y = refPos,
                operation = if(operation == 'S') operation else "M" # operation
            )
            qryPos <<- qryPos + (if(operation == 'S') 0 else length)
            
        # M or m operations present in reference and query
        } else {
            dt <- data.table(
                x = qryPos:(qryPos + length - 1),                
                y = refPos:(refPos + length - 1),
                operation = "M" # operation
            )
            refPos <<- refPos + length
            qryPos <<- qryPos + length
        }
        dt
    }))
    notS <- dt$operation != "S"
    if(strand == "-"){
        nNotS <- sum(notS)
        dt$y[notS] <- dt$y[notS][nNotS:1]
    }
    list(
        qryPos = qryPos,
        refPos = refPos,
        dt = dt,
        x = mean(dt$x[notS]),
        y = mean(dt$y[notS])
    )
}

# plot all alignment segments from a single read (could be 1 or 2 reads per moleculeType)
plotReadAlns <- function(pts, mmd, readN, chrom, midx){
    qryPos <- 1
    is <- which(mmd$alnReadNs == readN)
    types    <- mmd$types[is]
    chroms   <- mmd$chroms[is]
    strands  <- mmd$strands[is]
    poss     <- mmd$poss[is]
    cigars   <- mmd$cigars[is]
    mapQs    <- mmd$mapQs[is]
    baseQuals<- mmd$baseQuals[is]
    is <- which(mmd$jxnReadNs == readN)
    sizes    <- mmd$sizes[is]
    insSizes <- mmd$insSizes[is]
    for(i in 1:length(types)){
        refPos <- getRefPos(readN, strands[i], poss[i], cigars[i])
        dp <- cigarDotPlot(cigars[i], qryPos, refPos, strands[i])
        text(
            dp$x, dp$y, 
            paste(
                paste("QUAL", baseQuals[i], sep = " "), 
                paste("MAPQ", mapQs[i], sep = " "),
                sep = ", "
            ),
            pos = if(dp$x > midx) 2 else 4, 
            offset = 1.5
        ) 
        pts <- rbind(pts, data.table(
            x = dp$dt$x,
            y = dp$dt$y,
            col = unlist(dotPlotColors[dp$dt$operation]),
            operation = dp$dt$operation
        ))
        if(isTruthy(types[i]) && types[i] %in% c("D","I")){
            if(!is.na(sizes[i]) && sizes[i] > 0 && sizes[i] < midx * 2) pts <- rbind(pts, data.table(
                x = rep(dp$qryPos, sizes[i]),
                y = dp$refPos:(dp$refPos + sizes[i] - 1),
                col = if(chroms[i] == chrom) dotPlotColors$D else NA,
                operation = "D"
            ))     
            if(!is.na(insSizes[i]) && insSizes[i] > 0) pts <- rbind(pts, data.table(
                x = dp$qryPos:(dp$qryPos + insSizes[i] - 1),
                y = rep(dp$refPos + sizes[i], insSizes[i]),
                col = if(chroms[i] == chrom) dotPlotColors$I else NA,
                operation = "I"
            ))   
            if(!is.na(insSizes[i]) && insSizes[i] < 0) pts <- rbind(pts, data.table(
                x = (dp$qryPos + insSizes[i] + 1):dp$qryPos,
                y = rep(dp$refPos + sizes[i], -insSizes[i]),
                col = if(chroms[i] == chrom) dotPlotColors$H else NA,
                operation = "H"
            )) 
        }
        qryPos <- dp$qryPos + if(!is.na(insSizes[i])) insSizes[i] else 0  
    }
    pts     
}

# plot a single chromosome segment
renderAlignmentPlot <- function(xmax, mmd, i, ampliconic){
    ylim <- if(ampliconic) c(mmd$amplicon$pos1, mmd$amplicon$pos2) else {
        range(mapply(function(ampliconic, readN, strand, pos, cigar){
            if(!ampliconic) return(rep(as.integer(NA), 2))
            p <- getRefPos(readN, strand, pos, cigar)
            l <- getNRefBases(cigar)       
            c(p, p + l - 1)     
        }, mmd$ampliconic, mmd$alnReadNs[i], mmd$strand[i], mmd$pos[i], mmd$cigar[i]), na.rm = TRUE)
    }
    chrom <- if(ampliconic) mmd$amplicon$chrom1 else mmd$chroms[i][1]
    par(mar = c(if(ampliconic) 4.1 else 0.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = c(1, xmax),
        ylim = ylim,
        xlab = if(ampliconic) "Position in Read" else "",
        ylab = chrom,
        xaxt = if(ampliconic) "s" else "n"
    )
    pts <- data.table(x = integer(), y = integer(), col = character(), operation = character())
    midx <- xmax / 2 

    pts <- plotReadAlns(pts, mmd, 1, chrom, midx)    
    if(!mmd$isMerged) pts <- plotReadAlns(pts, mmd, 2, chrom, midx)
    pts <- pts[order(unlist(dotStackOrder[pts$operation]))]
    abline(v = midx, col = CONSTANTS$plotlyColors$black) 
    abline(v = c(midx + 50 * c(1:10, -(1:10))), col = CONSTANTS$plotlyColors$grey) 
    midy <- mean(ylim)
    abline(h = midy, col = CONSTANTS$plotlyColors$black)
    abline(h = c(midy + 50 * c(1:10, -(1:10))), col = CONSTANTS$plotlyColors$grey)
    points(pts$x, pts$y, col = pts$col, pch = 19, cex = 0.5)
}

# plot the read(s) QUAL profile(s)
renderReadQualPlot <- function(xmax, mmd, mt){
    par(mar = c(0.1, 4.1, 0.1, 0.1), cex = 1)
    plot(
        NA, NA,
        xlim = c(1, xmax),
        ylim = c(0, 42),
        ylab = "QUAL",
        xaxt = "n"
    )
    midx <- xmax / 2
    abline(v = midx, col = CONSTANTS$plotlyColors$black) 
    abline(v = c(midx + 50 * c(1:10, -(1:10))), col = CONSTANTS$plotlyColors$grey) 
    lines(
        x = 1:nchar(mt$seq1),
        y = as.integer(sapply(strsplit(mt$qual1, "")[[1]], charToRaw)) - 33,
        col = CONSTANTS$plotlyColors$blue
    )
    if(!mmd$isMerged) lines(
        x = 1:nchar(mt$seq2),
        y = as.integer(sapply(strsplit(mt$qual2, "")[[1]], charToRaw)) - 33,
        col = CONSTANTS$plotlyColors$orange
    )    
}

# render the composite QC plot; this is the main function call, cascading upwards
moleculeQcPlot <- function(moleculeMetadata) {
    mmd <- moleculeMetadata()
    req(mmd)
    mt <- mmd$moleculeType
    xmax <- max(mmd$ampliconSize, mmd$tLen)

    # prepare for out-of-amplicon rows
    if(any(!mmd$ampliconic)){
        nOutSegments <- sum(!mmd$ampliconic)
        nPlotRows <- nOutSegments + 2
        heights <- c(0.15, rep(0.15, nOutSegments), 0.85 - 0.15 * nOutSegments)
    } else {
        nPlotRows <- 2
        heights <- c(0.15, 0.85)
    }
    layout(matrix(1:nPlotRows, ncol = 1), heights = heights) 

    # plot QUAL
    renderReadQualPlot(xmax, mmd, mt)

    # plot out-of-amplicon alignments; at present, only handle one such segment per chromosome well
    for(i in which(!mmd$ampliconic)){
        renderAlignmentPlot(xmax, mmd, i, FALSE) # these could be from different regions!
    }

    # plot ampliconic alignments
    renderAlignmentPlot(xmax, mmd, which(mmd$ampliconic), TRUE)
}
