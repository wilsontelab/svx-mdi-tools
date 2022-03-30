
#----------------------------------------------------------------------
# show_junction.R creates "zoomed" representations of a single SV junction
#----------------------------------------------------------------------

# some common plot values
jxn.cex.text = 1 # 1.5
jxn.cex.legend = 1 # 1.25
jxn.cex.point = 1.1 # 1.25
maxJxnDist <- 1000

# return information on the SV we are zooming in on
zoomInfo <- reactive({
    reportProgress('getZoomInfo')
    if(input$sample_1 == '-' | is.null(input$svTableOutput_rows_selected)) return(NULL)
    rowI <- as.integer(input$svTableOutput_rows_selected)
    smp  <- svTable[rowI,smp]
    svId <- svTable[rowI,svId]
    nodes <- getSvMoleculesFromList(smp, svId)
    nodes[,c('chrom','side','pos') := unpackNodeNames(NODE)]
    list(
        smp  = smp,
        svId = svId,
        sv = sampleData[[smp]]$svTable[SV_ID == svId],
        nodes = nodes[order(side,chrom,pos)] # thus, nodes are in canonical order
    )
})
isSequencedJunction <- function(sv) sv[,!is.na(JXN_SEQ)] # either sequenced or reconstructed
isReconstructedJunction <- function(sv) sv[,!is.na(JXN_SEQ) & !is.na(MERGE_LEN)]

#----------------------------------------------------------------------
# plot of stacked read sequences with base-color alignment visualization
#----------------------------------------------------------------------
wsh <- 0.7
zoomBaseColors <- list( # generally follow IGV base color conventions
    A = c(0,  1,  0), # green 
    a = c(wsh,1,  wsh), 
    C = c(0,  0,  1), # blue
    c = c(wsh,wsh,1),
    G = c(0.82,0.43,0),
    g = c(0.94,0.83,wsh),
    T = c(1,  0,  0), # red
    t = c(1,  wsh,wsh),
    N = c(0,  0,  0), # black
    n = c(wsh,wsh,wsh),
    x = c(wsh,wsh,wsh), # clip mask
    '-' = c(0.9,0.9,0.9), # deleted/missing = light grey
    '+' = c(1,0,1), # insertion = purple
    ' ' = c(1,1,1) # white-space
)
clipBases <- c('a','c','g','t','n')
zoomBaseR <- sapply(zoomBaseColors, function(v) v[1])
zoomBaseG <- sapply(zoomBaseColors, function(v) v[2])
zoomBaseB <- sapply(zoomBaseColors, function(v) v[3])
#-------------------------------------------------------------------------------------
unpackNodeName <- function(nodeName){ # extract chrom, side, pos from node sort keys, e.g. chr1:R:32538
    x <- as.list(strsplit(nodeName, '\\W')[[1]])
    x[[3]] <- as.integer(x[[3]])
    names(x) <- c('chrom','side','pos')
    x
}
unpackNodeNames <- function(nodeNames){ # extract chrom, side, pos from node sort keys, e.g. chr1:R:32538
    x <- as.data.frame(
        matrix(unlist(strsplit(nodeNames, '\\W')), ncol=3, byrow=TRUE),
        stringsAsFactors=FALSE
    )
    x[[3]] <- as.integer(x[[3]])
    x
}
consensusBase <- function(x){ # get the most frequent sequence value at a reference position from a character vector
    x <- na.omit(x)           # not necessarily a single base, e.g. could be "GT" for an insertion relative to reference
    if(length(x) == 0) " "
    else if(all(x == x[1])) x[1]
    else {
        x <- sort(table(x),decreasing=TRUE)
        if(x[1] == x[2]) "N" else names(x[1])
    }
}
jxnImgData <- reactiveVal() # for cascading data from createJunctionZoom to createJunctionConsensus
createJunctionConsensus <- function(d, zoom, png, clipMode){
    makeTrackImg <- function(bases){ # helper to create an image row for reference or consensus
        bases <- t(sapply(bases, function(base) zoomBaseColors[[base]]))
        array(bases, dim=c(nrow(bases),1,3))  
    }
    isInsPos <- if(zoom$sv[1,MICROHOM_LEN] >= 0) FALSE else {
        nCol <- ncol(d$baseMap) # create a matrix for protecting insertion bases from being masked out of consensus
        t(apply(d$baseMap, 1, function(v) { 
            isIns <- all(na.omit(v) %in% clipBases)
            rep(isIns, nCol)
        }))
    }
    
    # make the consensus track
    consensus <- matrix(ifelse(d$baseMap %in% clipBases & !isInsPos, NA, d$baseMap), nrow=nrow(d$baseMap))
    consensus <- apply(consensus, 1, consensusBase)
    consensus <- makeTrackImg(consensus)
    
    # make the reference track matching left side of junction ...
    isRc1 <- zoom$sv$SIDE_1 == "R"
    if(isRc1) zoom$sv$GEN_REF_1 <- rc(zoom$sv$GEN_REF_1)
    ref1All <- strsplit(zoom$sv$GEN_REF_1, '')[[1]] # the sequence of all reference bases within the sequenced window
    ref1Masked <- ifelse(1:length(ref1All) <= (zoom$sv[1,FAIDX_PADDING]+1), ref1All, " ") # now with " " for positions past the junctions
    ref1 <- makeTrackImg(ref1Masked[d$usedPos])
    
    # ... and the right side of junction
    isRc2 <- zoom$sv$SIDE_2 == "L"
    if(isRc2) zoom$sv$GEN_REF_2 <- rc(zoom$sv$GEN_REF_2)
    ref2All <- strsplit(zoom$sv$GEN_REF_2, '')[[1]]
    ref2Masked <- ifelse(1:length(ref2All) > (length(ref2All) - (zoom$sv[1,FAIDX_PADDING]+1)), ref2All, " ")
    ref2 <- makeTrackImg(ref2Masked[d$usedPos + zoom$sv[1,MICROHOM_LEN] - 1])
    
    # put the three tracks together, save the image are return its dimensions
    imgMap <- abind(ref1, consensus, ref2, along=2)
    imgMap[is.na(imgMap)] <- 1
    suppressWarnings(save.image(as.cimg(imgMap), png))
    dim(imgMap)
}
createJunctionZoom <- function(zoom, png, clipMode){

    # at first, only plot the nodes at the junction (can add proper gaps later)
    zoom$nodes <- zoom$nodes[IS_JUNCTION_NODE==TRUE] 
    
    # initialize image and dimensions
    zoomMapPadding <- zoom$sv[1,FAIDX_PADDING] # set width used when preparing the output (wider than final output)
    leftRefI <- zoomMapPadding + 1
    microhomologyLength <- zoom$sv[1,MICROHOM_LEN]
    zoomMapWidth <- 2 * zoomMapPadding + 2 - microhomologyLength
    nAln <- zoom$nodes[,.N] # the number of rows needed for graphical display
    vPadding <- 3
    yOffset <- 0
    baseMap <- matrix(NA, nrow=zoomMapWidth, ncol=nAln) # for consensus making
    imgMap  <- array(NA, dim=c(zoomMapWidth, nAln+vPadding*2, 3)) # for graphical display

    # set verticals to aid in junction localization
    if(zoom$sv[,MICROHOM_LEN] == 0){
        imgMap[leftRefI,,] <- 0 # demarcate stem bases 1-base INSIDE any microhomology/insertion            
        imgMap[leftRefI+1,,] <- 0                
    } else if(zoom$sv[,MICROHOM_LEN] > 0){ # micromology = yellow
        left  <- leftRefI - zoom$sv[,MICROHOM_LEN] + 1
        right <- leftRefI
        imgMap[left-1,,] <- 0 # demarcate stem bases 1-base INSIDE any microhomology/insertion            
        imgMap[right+1,,] <- 0  
        is <- left:right
        imgMap[is,,1] <- 0.9
        imgMap[is,,2] <- 0.9
        imgMap[is,,3] <- 0
    } else if(zoom$sv[,MICROHOM_LEN] < 0){
        left  <- leftRefI + 1
        right <- leftRefI - zoom$sv[,MICROHOM_LEN]
        imgMap[left-1,,] <- 0 # demarcate stem bases 1-base INSIDE any microhomology/insertion            
        imgMap[right+1,,] <- 0  
        is <- left:right
        imgMap[is,,1] <- 1
        imgMap[is,,2] <- 0
        imgMap[is,,3] <- 1
    }    
    
    # add all alignments to maps
    leftNodeN <- zoom$nodes[IS_REF_NODE==TRUE][1,NODE_N]
    for (nodeN in if(leftNodeN==1) 1:2 else 2:1){
        isLeftNode <- nodeN == leftNodeN
        baseMapNodeI <- if(isLeftNode) 1 else 2
        sideNodes <- zoom$nodes[NODE_N==nodeN][order(if(isLeftNode) CLIP_LEN else -CLIP_LEN)]
        nSideAln <- sideNodes[,.N]
        refNode <- sideNodes[IS_REF_NODE==TRUE,unpackNodeName(NODE)]
        for(i in 1:nSideAln){   
            node <- sideNodes[i,unpackNodeName(NODE)]
            parsedAln <- cigarToRefAln(sideNodes[i,CIGAR], sideNodes[i,SEQ])            
                 if(clipMode == "grey") parsedAln$seq <- ifelse(parsedAln$seq %in% clipBases, 'x', parsedAln$seq)
            else if(clipMode == "none") parsedAln$seq <- ifelse(parsedAln$seq %in% clipBases, ' ', parsedAln$seq)
            nodePosDelta <- if(node$side == "L"){
                refNode$pos - node$pos
            } else {
                node$pos - refNode$pos
            }
            mapPos <- if(isLeftNode){
                zoomMapPadding + 1 - nodePosDelta - parsedAln$lengthAlignedOut - parsedAln$leftClip + 1
            } else {
                zoomMapPadding + 1 - zoom$sv[,MICROHOM_LEN] + 1 + nodePosDelta - parsedAln$leftClip
            }
            j <- mapPos:(mapPos + parsedAln$lengthOut - 1)
            col <- t(sapply(parsedAln$seq, function(base) zoomBaseColors[[base]]))
            imgMap[j,i+yOffset+vPadding,] <- col
            baseMap[j,i+yOffset] <- parsedAln$seq
        }
        yOffset <- nSideAln
    }
    
    # assemble the final image
    usedPos <- rle(apply(baseMap[,], 1, function(v) any(!is.na(v))))$lengths
    usedPos <- (usedPos[1] + 1):(zoomMapWidth - usedPos[length(usedPos)]) # row indices = base positions in SV allele
    baseMap <- baseMap[usedPos,]
    imgMap <- imgMap[usedPos,,]
    imgMap[is.na(imgMap)] <- 1
    suppressWarnings(save.image(as.cimg(imgMap), png))
    
    # pass data on to consensus maker
    jxnImgData(list(usedPos=usedPos, baseMap=baseMap))
    
    # return the dimension of the read map
    dim(imgMap)
}

#----------------------------------------------------------------------
# text-based sequence view of a single selected junction
#----------------------------------------------------------------------
makeJunctionHtml <- function(zoom){
    reportProgress('makeJunctionHtml')
    
    # determine the requested alignment type
    if(input$jxnDisplayType == "Genome"){
        offset <- zoom$idx$ALIGNS_OFFSET
        length <- zoom$idx$ALIGNS_LENGTH
    } else if(input$jxnDisplayType == "Genome+") {
        offset <- zoom$idx$ALIGNS_PLUS_OFFSET
        length <- zoom$idx$ALIGNS_PLUS_LENGTH
    } else if(input$jxnDisplayType == "Reads") {
        offset <- zoom$idx$READS_OFFSET
        length <- zoom$idx$READS_LENGTH        
    }
    
    # get the alignment data
    data <- getSvIndexed(zoom$smp, zoom$svId, 'alignments', offset, length)
    
    # parse into readable rows
    nLines <- length(data)
    nChar <- nchar(data[1])
    charPerLine <- as.numeric(input$charPerLine)
    nChunks <- ceiling(nChar / charPerLine)
    nCharLastLine <- nChar %% charPerLine
    if(nCharLastLine==0) nCharLastLine <- charPerLine
    output <- character()
    blueLines <- list( # junction lines to denote with highlight coloring
        Genome = 3:4,
        'Genome+' = 3:4,
        Reads = 1  
    )
    for(i in 1:nChunks){
        if(i > 1) output <- c(output, "<br>")        
        start <- 1 + (i - 1) * charPerLine
        for(j in 1:nLines){
            str  <- substr(data[j], start, start + charPerLine - 1)
            str  <- gsub(' ', '&nbsp;', str) 
            str_ <- gsub('~', '', str) 
            if(nchar(str_) > 0){
                col <- if(j %in% blueLines[[input$jxnDisplayType]]) 'blue' else 'black'
                output <- c(
                    output,
                    paste(paste('<span style="color:', col, ';">', sep=""),
                          str,
                          '</span>', sep="")
                )                
            }            
        }
    }
    output <- gsub('~', '&nbsp;', output)
    
    # app/prepend the outermost chromosomal positions
    firstInPair <- 64
    svx <- if(bitwAnd(zoom$sv$FLAG1, firstInPair)) {
        zoom$sv[,.(RNAME1,PROX_OUT_POS1,RNAME2,PROX_OUT_POS2)]
    } else {
        zoom$sv[,.(RNAME2,PROX_OUT_POS2,RNAME1,PROX_OUT_POS1)]
    }
    output <- c(
        "<br>",
        paste(svx[[1]], format(svx[[2]], big.mark=","), sep=":"),
        '&darr;',
        output,
        paste(paste(rep("&nbsp;", nCharLastLine-1),collapse=""),'&uarr;',sep=""),
        paste(svx[[3]], format(svx[[4]], big.mark=","), sep=":")
    )
    
    # return the hopefully readable sequence block
    paste(output, collapse="<br>")
}

#----------------------------------------------------------------------
# plot size by placement of molecules supporting a single selected junction
#----------------------------------------------------------------------
getSvEvidence <- function(zoom, fileType, offset, length){ # get the source molecules data for an SV
    reportProgress('getSvEvidence')   
    data <- getSvIndexed(zoom$smp, zoom$svId,
                         fileType, zoom$idx[[offset]], zoom$idx[[length]])
    if(is.null(data)) return() 
    data <- data.frame(
        matrix(unlist(strsplit(data, "\t")), length(data), length(find[[fileType]]), byrow=TRUE),
        stringsAsFactors=FALSE
    )
    colnames(data) <- names(find[[fileType]])
    rownames(data) <- 1:nrow(data)
    data
}
getJxnDist <- function(mol, zoom, jxnSide, molSide=jxnSide){
    jxnPos <- zoom$sv[[paste('PROX_JXN_POS', jxnSide, sep="")]]
    molOutPos <- as.integer(mol[paste('MOL_OUT_POS', molSide, sep="")]) # non-zero for split reads
    if(molOutPos == 0) molOutPos <- as.integer(mol[paste('OUT_POS', molSide, sep="")])    
    dist <- max(molOutPos,jxnPos) - min(molOutPos,jxnPos) + 1 - zoom$sv$MICROHOM_LEN / 2
    if(dist <= maxJxnDist) dist else NA  # may get wacky values for molecules with multiple junctions
}
getClipDist <- function(mol, zoom, jxnSide, molSide=jxnSide){
    jxnPos <- zoom$sv[[paste('PROX_JXN_POS', jxnSide, sep="")]]
    molOutPos <- as.integer(mol[paste('OUT_POS', molSide, sep="")])    
    max(molOutPos,jxnPos) - min(molOutPos,jxnPos) + 1 - zoom$sv$MICROHOM_LEN / 2
}
addClipsToPlot <- function(zoom){

    # load clip data (all columns still as strings here)
    clips <- getSvEvidence(zoom, 'clips', 'CLIP_OFFSET', 'CLIP_LENGTH')
    if(is.null(clips)) return()
    
    # convert required columns to integers
    for (col in c('OUT_POS1','OUT_POS2','IS_DUPLEX','OUT_CLIP1')) {
        clips[[col]] <- as.integer(clips[[col]])
    }
    
    # calculate plot values
    jxnSide <- ifelse(clips$OUT_POS1 == zoom$sv$PROX_JXN_POS1, 1, 2)
    dist1 <- ifelse(jxnSide == 1, apply(clips, 1, getClipDist, zoom, 1, 2), clips$OUT_CLIP1)    
    dist2 <- ifelse(jxnSide == 1, clips$OUT_CLIP1, apply(clips, 1, getClipDist, zoom, 2, 2))
    size <- dist1 + dist2
    rpos <- dist1 - dist2 
    
    # add clip data points
    points(
        rpos, size,
        pch=4, # X = clip
        col=ifelse(clips$IS_DUPLEX==1, plotColors$duplex$duplex, plotColors$duplex$single)
    )
}
plotMolSizePos <- function(zoom){

    # load molecule data (all columns still as strings here)
    mol <- getSvEvidence(zoom, 'molecules', 'MOL_OFFSET', 'MOL_LENGTH')
    
    # calculate the distance from molecule outer ends to junctions (adjusted for microhomology length)
    dist1 <- apply(mol, 1, getJxnDist, zoom, 1)
    dist2 <- apply(mol, 1, getJxnDist, zoom, 2)
    
    # calculate the inferred size of molecule and the junction position within it
    size <- dist1 + dist2
    rpos <- dist1 - dist2 # so deletion maps left to right across SV fragment
    
    # set plot symbols and colors
    mol$pch <- ifelse(
        mol$JXN_CLASS==jxnClasses$GAP,
        1,           # open circle = unsequenced junction in unmerged molecule gap   
        ifelse(mol$IS_MERGED=='1',
               16,   # close circle = sequenced junction in merged molecule
               15)) # sequenced junction in one read of an unmerged molecule
    mol$col <- ifelse(mol$IS_DUPLEX==1, plotColors$duplex$duplex, plotColors$duplex$single)
    
    # initialize the plot
    sd <- sampleData[[zoom$smp]]
    cex_ <- as.numeric(input$pointCex)
    plot(
        0, 0, typ="n",
        xlim=c(-sd$p95,sd$p95), ylim=c(0,sd$p99),
        xlab="Relative position in molecule",ylab="Mol. size (bp)",
        cex.axis=jxn.cex.text * cex_,
        cex.lab=jxn.cex.text * cex_
    )    
    
    # add line rules
    abline(h=projectInfo$READ_LENGTH, col="black")
    abline(h=c(sd$p5,sd$p50,sd$p95), col=c("darkgrey","black","darkgrey"))
    abline(v=0, col="darkgrey")
    abline(0,  1, col="darkgrey")
    abline(0, -1, col="darkgrey")
    
    # add legends
    legend("bottomleft",
           c('seq by read','seq merged','in gap','clipped'),
           pch=c(15,16,1,4),
           cex=jxn.cex.legend * cex_)
    legend("bottomright",
           names(plotColors$duplex),
           col=unlist(plotColors$duplex),
           pch=16,
           cex=jxn.cex.legend * cex_)
    
    # add one data point for each source molecule (GAP or SPLIT)
    # symbols and colors denote the nature of the molecule evidence
    points(rpos, size, pch=mol$pch, col=mol$col, cex=jxn.cex.point * cex_)
    
    # since find.molecules not reporting CLIPS, add them here from tabixed extract file
    addClipsToPlot(zoom)
}


# THIS CODE CREATES A PLOT TO HELP UNDERSTAND THE PATTERN EXPECTED
#minMappable <- 30
#minMappable_2x <- 2 * minMappable
#readLen <- 150
#minMergeOverlap <- 10
#maxMergedSize <- 2 * readLen - minMergeOverlap
#maxTLen <- 500
#
#nMol <- 2500
#molSize <- sample(1:maxTLen, nMol, replace=TRUE)
#jxnPos <- round(pmin(pmax(runif(nMol, 0, 1) * molSize, 1), molSize),0)
#
#x <- jxnPos - molSize / 2
#
#shortSide <- pmin(jxnPos, molSize - jxnPos)
#longSide <- pmax(jxnPos, molSize - jxnPos)
#isMappable <- longSide >= minMappable
#isMerged <- molSize < maxMergedSize
#isOuterClip <- shortSide < minMappable
#inGap <- !isMerged & shortSide > readLen
#isInnerClip <- !isMerged & !inGap & shortSide > (readLen - minMappable)
#
#col <- ifelse(!isMappable, NA,
#    ifelse(isOuterClip, "black",
#        ifelse(isInnerClip, "orange3",
#            ifelse(isMerged, "green4",
#                "blue"
#            )
#        )
#    )
#)
#pch <- ifelse(inGap, 1, 16)
#plot(jxnPos,molSize-jxnPos,xaxs="i",yaxs="i", pch=pch, col=col)
#abline(0,1,lwd=2,col="grey")
#abline(maxMergedSize,-1,lwd=2,col="grey")
#h <- c(minMappable,readLen)
#abline(h=h,v=h)
#
#
#plot(x,molSize,xaxs="i",yaxs="i", pch=pch, col=col)
#abline(0,2)
#abline(0,-2)
#abline(h=c(minMappable,minMappable_2x,maxMergedSize))
#abline(v=0)





    #clipFile <- getExtractFile(zoom$smp, 'sortedClips.bgz')
    #addClips <- function(jxnSide){
    #    
    #    # get sv values
    #    chrom    <- zoom$sv[[paste('RNAME', jxnSide, sep="")]]
    #    jxnPos   <- zoom$sv[[paste('PROX_JXN_POS', jxnSide, sep="")]]
    #    
    #    # pull clips from tabixed file
    #    region <- paste(chrom, ":", jxnPos, "-", jxnPos, sep="")
    #    pipe <- pipe(paste('tabix', clipFile, region))
    #    clips <- read.table(
    #        pipe,
    #        sep="\t",
    #        header=FALSE,
    #        stringsAsFactors=FALSE,
    #        col.names = names(extract$clips),
    #        colClasses = extract$clips
    #    )
    #    
    #    # restrict to those that match the expected junction direction        
    #    jxnSideX <- zoom$sv[[paste('INN_SIDE', jxnSide, sep="")]]
    #    clips <- clips[clips$OUT_SIDE1==jxnSideX,]
    #    
    #    # calculate plot values
    #    if(jxnSide == 1){
    #        dist1 <- apply(clips, 1, getJxnDist, zoom, jxnSide, 2)
    #        dist2 <- clips$OUT_CLIP1
    #    } else {
    #        dist1 <- clips$OUT_CLIP1
    #        dist2 <- apply(clips, 1, getJxnDist, zoom, jxnSide, 2)  
    #    }
    #    size <- dist1 + dist2
    #    rpos <- (dist2 - dist1) 
    #    
    #    # add clip data points
    #    points(
    #        rpos, size,
    #        pch=4, # X = clip
    #        col=ifelse(clips$IS_DUPLEX==1, plotColors$duplex$duplex, plotColors$duplex$single)
    #    )
    #}
    #addClips(1)
    #addClips(2)

# data.table order (radix)
    #   a     x
    #1: 1  chr1
    #2: 1 chr10
    #3: 1 chr11
    #4: 1 chr12
    #5: 1  chr2
    #6: 1 chr20
    #7: 1  chrX
# linux sort
    #(base) [wilsonte@wilsonte-n1 svCapture]$ echo $x | sed 's/ /\n/g' | sort -k1,1
    #chr1
    #chr10
    #chr11
    #chr12
    #chr2
    #chr20
    #chrX
# perl sort
    # perl -e '@x=split(" ","chr1 chr2 chr10 chr11 chr12 chr20 chrX"); print join("\n",sort {$a cmp $b} @x),"\n" '
    #chr1
    #chr10
    #chr11
    #chr12
    #chr2
    #chr20
    #chrX

# data.table order (radix)
    #> d[order(x)]
    #   a       x
    #1: 1 chr10:L
    #2: 1 chr11:L
    #3: 1 chr12:L
    #4: 1  chr1:L
    #5: 1 chr20:L
    #6: 1  chr2:L
    #7: 1  chrX:L
# linux sort
    #(base) [wilsonte@wilsonte-n1 svCapture]$ echo $y | sed 's/ /\n/g' | sort -k1,1
    #chr10:L
    #chr11:L
    #chr12:L
    #chr1:L
    #chr20:L
    #chr2:L
    #chrX:L
# perl sort
    #(base) [wilsonte@wilsonte-n1 svCapture]$ perl -e '@x=split(" ","chr1:L chr2:L chr10:L chr11:L chr12:L chr20:L chrX:L"); print join("\n",sort {$a cmp $b} @x),"\n" '
    #chr10:L
    #chr11:L
    #chr12:L
    #chr1:L
    #chr20:L
    #chr2:L
    #chrX:L
