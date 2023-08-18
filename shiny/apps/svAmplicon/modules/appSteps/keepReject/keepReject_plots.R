#----------------------------------------------------------------------
# construct a composite paired quality plot, i.e. MAPQ vs. QUAL
#----------------------------------------------------------------------
pairedQualityPlotServer <- function(moleculeTypes) mdiInteractivePlotServer(
    "pairedQualityPlot",       
    click = TRUE,
    contents = reactive({ 
        mts <- moleculeTypes()        
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
        refPos <- min(dt$y[notS])
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
plotReadAlns <- function(pts, mmd, readN, midx){
    qryPos <- 1
    alnIs <- 1:length(mmd$chroms)
    nodeIs <- sapply(alnIs * 2, function(j) (j - 1):j)
    for(alnI in alnIs){
        chr <- mmd$chroms[alnI]
        nodeJ <- alnI * 2
        refPos <- min(mmd$poss[(nodeJ - 1):nodeJ])
        dp <- cigarDotPlot(mmd$cigars[alnI], qryPos, refPos, mmd$strands[alnI]) 
        pts <- rbind(pts, data.table(
            x = dp$dt$x,
            y = dp$dt$y,
            col = unlist(dotPlotColors[dp$dt$operation]),
            operation = dp$dt$operation,
            chrom = chr,
            alnI = alnI
        ))
        eT <- mmd$edgeTypes[alnI]
        if(isTruthy(eT)){
            chr <- mmd$chroms[alnI]
            eS  <- mmd$eventSizes[alnI]
            iS  <- mmd$insertSizes[alnI]
            if(!is.na(eS) && eS > 0 && eS < midx * 2 && eT %in% c("D","U","V")) pts <- rbind(pts, data.table( # deletions / duplications / inversions
                x = rep(dp$qryPos, eS),
                y = switch( eT,
                    D = dp$refPos:(dp$refPos + eS - 1),
                    U = dp$refPos:(dp$refPos - eS + 1),
                    V = dp$refPos:(dp$refPos + eS - 1)
                ),
                col = dotPlotColors$D,
                operation = "D",
                chrom = chr,
                alnI = alnI
            ))     
            if(!is.na(iS) && iS > 0) pts <- rbind(pts, data.table( # junction insertion
                x = dp$qryPos:(dp$qryPos + iS - 1),
                y = switch(eT,
                    D = rep(dp$refPos + eS - 1, iS),
                    U = rep(dp$refPos - eS + 1, iS),
                    V = rep(dp$refPos + eS, iS),
                    rep(dp$refPos, iS)
                ),
                col = dotPlotColors$I,
                operation = "I",
                chrom = chr,
                alnI = alnI
            ))   
            if(!is.na(iS) && iS < 0) pts <- rbind(pts, data.table( # junction microhomology
                x = (dp$qryPos + iS + 1):dp$qryPos,
                y = switch(eT,
                    D = rep(dp$refPos + eS - 1, -iS),
                    U = rep(dp$refPos - eS + 1, -iS),
                    V = rep(dp$refPos + eS, -iS),
                    rep(dp$refPos, -iS)
                ),
                col = dotPlotColors$H,
                operation = "H",
                chrom = chr,
                alnI = alnI
            )) 
        }
        qryPos <- dp$qryPos
        if(!is.na(mmd$insertSizes[alnI])){
            qryPos <- qryPos + mmd$insertSizes[alnI]
        }
    }
    pts     
}

# plot a single chromosome segment
initializeAlignmentPlot <- function(xmax, mmd, alnI, ampliconic){
    ylim <- if(ampliconic) c(mmd$amplicon$pos1, mmd$amplicon$pos2) else { # the reference range on chrom for y axis
        j <- alnI * 2
        range(mmd$poss[(j - 1):j])
    }    
    par(mar = c(if(ampliconic) 4.1 else 0.1, 4.1, 0.1, 0.1), cex = 1)
    chrom_ <- if(ampliconic) mmd$amplicon$chrom1 else mmd$chroms[alnI]
    plot(
        NA, NA,
        xlim = c(1, xmax),
        ylim = ylim,
        xlab = if(ampliconic) "Position in Read" else "",
        ylab = chrom_,
        xaxt = if(ampliconic) "s" else "n"
    )
    if(mmd$isMerged){
        len <-mmd$tLen
        ovlp <- mmd$moleculeType$overlap
        flank <- (len - ovlp) / 2
        rect(max(0, flank), ylim[1], min(len, flank + ovlp), ylim[2], col = "grey90", lwd = NA) 
    } 
    midx <- xmax / 2     
    midy <- mean(ylim)    
    abline(v = midx, col = CONSTANTS$plotlyColors$black) 
    abline(v = c(midx + 50 * c(1:10, -(1:10))), col = CONSTANTS$plotlyColors$grey) 
    abline(h = midy, col = CONSTANTS$plotlyColors$black)
    abline(h = c(midy + 50 * c(1:10, -(1:10))), col = CONSTANTS$plotlyColors$grey)    
}
renderAlignmentPoints <- function(mmd, alnI_, ampliconic, pts){
    chrom_ <- if(ampliconic) mmd$amplicon$chrom1 else mmd$chroms[alnI_]
    pts <- pts[chrom == chrom_ & alnI == alnI_]
    pts <- pts[order(unlist(dotStackOrder[pts$operation]))]
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
        x = 1:nchar(mt$seq2) + mt$tLen1,
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

    # calculate plot points
    pts <- data.table(x = integer(), y = integer(), col = character(), operation = character(), 
                      chrom = character(), alnI = integer())
    midx <- xmax / 2 
    pts <- plotReadAlns(pts, mmd, 1, midx)    
    if(!mmd$isMerged) {  
        pts <- plotReadAlns(pts, mmd, 2, midx)
    }

    # plot out-of-amplicon alignments; at present, only handle one such segment per chromosome well
    for(alnI in which(!mmd$ampliconic)){
        initializeAlignmentPlot(xmax, mmd, alnI, FALSE)
        renderAlignmentPoints(mmd, alnI, FALSE, pts)
    }

    # plot ampliconic alignments
    alnIs <- which(mmd$ampliconic)
    initializeAlignmentPlot(xmax, mmd, alnIs[1], TRUE)
    for(alnI in alnIs){
        renderAlignmentPoints(mmd, alnI, TRUE, pts) # these could be from different regions!
    }
}
