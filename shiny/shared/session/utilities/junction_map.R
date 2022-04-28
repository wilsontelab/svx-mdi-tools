# ----------------------------------------------------------------------
# map of all base positions at and around a selected SV junction
# supports consensus making, evidence plotting, etc.
# ----------------------------------------------------------------------
wsh <- 0.7 # for washed-out colors
CONSTANTS$junctionBaseColors <- list( # generally follow IGV base color conventions
    A = c(0,    1,    0), # green 
    a = c(wsh,  1,    wsh), 
    C = c(0,    0,    1), # blue
    c = c(wsh,  wsh,  1),
    G = c(0.82, 0.43, 0),
    g = c(0.94, 0.83, wsh),
    T = c(1,    0,    0), # red
    t = c(1,    wsh,  wsh),
    N = c(0,    0,    0), # black
    n = c(wsh,  wsh,  wsh),
    x = c(wsh,  wsh,  wsh), # clip mask
    '-' = c(0.9, 0.9, 0.9), # deleted/missing = light grey
    '+' = c(1,   0,   1),   # insertion = purple
    ' ' = c(1,   1,   1)    # white-space
)
CONSTANTS$clipBases <- c('a', 'c', 'g', 't', 'n')

# convert an SV and its set of evidence molecules to a base map
getJunctionMap <- function(x, clipMode = "Faded Colors"){
    req(x)
    req(x$sv$JXN_BASES != "*")
    startSpinner(session, 'getJunctionMap') 

    # initialize map and dimensions
    refWidth <- nchar(x$sv$GEN_REF_1)
    faidxPadding <- (refWidth - 1) / 2
    microhomologyLength <- x$sv[, MICROHOM_LEN]    
    leftRefI  <- faidxPadding + 1
    rightRefI <- leftRefI + 1 - microhomologyLength
    mapWidth <- refWidth + 1 - microhomologyLength
    nAln <- x$mols[, .N + sum(!isOuterClip)]
    vPadding <- 3
    yOffset <- 1
    baseMap <- matrix(NA, nrow = mapWidth, ncol = nAln) # for consensus making
    imgMap  <- array(NA, dim = c(mapWidth, nAln + vPadding * 2, 3)) # for graphical display
    refMolAlns <- logical()

    # parse the genome reference lines
    GEN_REF_1 <- c(    strsplit(x$sv$GEN_REF_1, "")[[1]],  rep(" ", mapWidth))
    GEN_REF_2 <- c(rev(strsplit(x$sv$GEN_REF_2, "")[[1]]), rep(" ", mapWidth))
    GEN_REF_1 <- GEN_REF_1[1:mapWidth]
    GEN_REF_2 <- GEN_REF_2[1:mapWidth]
    GEN_REF_2 <- rev(GEN_REF_2)

    dprint(faidxPadding)
    dprint(nchar(x$sv$GEN_REF_1))
    dprint(nchar(x$sv$GEN_REF_2))
    dprint(x$sv$GEN_REF_1)
    dprint(x$sv$GEN_REF_2) 
    test <- strsplit(x$sv$GEN_REF_2, "")[[1]]
    dprint(paste(test[(faidxPadding + 1):(faidxPadding + 100)], collapse = ""))


    # set verticals as a visual aid to junction localization
    if(microhomologyLength == 0){
        imgMap[leftRefI, , ] <- 0 # demarcate two stem bases around a blunt joint   
        imgMap[leftRefI + 1, , ] <- 0
    } else if(microhomologyLength > 0){ # micromology = yellow
        left  <- leftRefI - microhomologyLength + 1
        right <- leftRefI
        imgMap[left - 1, , ] <- 0     
        imgMap[right + 1, , ] <- 0  
        is <- left:right
        imgMap[is, , 1] <- 0.9
        imgMap[is, , 2] <- 0.9
        imgMap[is, , 3] <- 0
    } else if(microhomologyLength < 0){
        left  <- leftRefI + 1
        right <- leftRefI - microhomologyLength
        imgMap[left - 1, , ] <- 0      
        imgMap[right + 1, , ] <- 0  
        is <- left:right
        imgMap[is, , 1] <- 1
        imgMap[is, , 2] <- 0
        imgMap[is, , 3] <- 1
    }    
    
    # add all alignments to maps
    processNode <- function(nodeN, isRef, jxnSide, jxnPos, CIGAR, SEQ, pos){
        if(CIGAR == "*") return(NULL) # the missing node of an outer clip
        refMolAlns <<- c(refMolAlns, isRef == 1)
        aln <- cigarToRefAln(CIGAR, SEQ)           
        nodePosDelta <- if(jxnSide == "L") jxnPos - pos else pos - jxnPos
        mapPos <- if(nodeN == 1){
            leftRefI - nodePosDelta - aln$lengthAlignedOut - aln$leftClip + 1
        } else {
            leftRefI + 1 - microhomologyLength + nodePosDelta - aln$leftClip
        }
        j <- mapPos:(mapPos + aln$lengthOut - 1)
        baseMap[j, yOffset] <<- aln$seq        
             if(clipMode == "Greyed Out") aln$seq[aln$seq %in% CONSTANTS$clipBases] <- "x"
        else if(clipMode == "Hidden")     aln$seq[aln$seq %in% CONSTANTS$clipBases] <- " "
        col <- t(sapply(aln$seq, function(base) CONSTANTS$junctionBaseColors[[base]]))        
        imgMap[ j, yOffset + vPadding, ] <<- col        
        yOffset <<- yOffset + 1 
    }
    x$mols[
        order( CLIP_LEN_1, POS_1), 
        mapply(processNode, 1, IS_REFERENCE, x$sv$SIDE_1, x$sv$POS_1, CIGAR_1, SEQ_1, pos1)
    ]
    x$mols[
        order(-CLIP_LEN_2, POS_2 + nchar(SEQ_2)), 
        mapply(processNode, 2, IS_REFERENCE, x$sv$SIDE_2, x$sv$POS_2, CIGAR_2, SEQ_2, pos2)
    ]

    # assemble the final maps
    usedPos <- rle(apply(baseMap[, ], 1, function(v) any(!is.na(v))))$lengths
    usedPos <- (usedPos[1] + 1):(mapWidth - usedPos[length(usedPos)]) # row indices = base positions in SV allele
    imgMap[is.na(imgMap)] <- 1

    # finish and return the results
    stopSpinner(session, 'getJunctionMap') 
    list(
        sv      = x$sv,
        mols    = x$mols,  
        text    = baseMap, 
        image   = imgMap,
        refMolAlns = refMolAlns,
        leftRefI = leftRefI,
        rightRefI = rightRefI,
        usedPos = usedPos,
        GEN_REF_1 = GEN_REF_1,
        GEN_REF_2 = GEN_REF_2
    )
}

# create an alignment betweeen junction molecules, their consensus and the genome reference
getJunctionAlignment <- function(map, charPerLine = 100, mode = "Evidence Consensus"){
    req(map)    

    # consensus mode
    if(mode == "Evidence Consensus"){
        consensus <- apply(map$text, 1, function(x){
            x <- x[!is.na(x)]
            if(length(x) == 0) return("~")
            xx <- x[!(x %in% CONSTANTS$clipBases)] # insertions are also lower case
            if(length(xx) == 0) xx <- x
            agg <- aggregate(xx, list(xx), length)
            agg[which.max(agg[[2]]), 1]
        })
        match1 <- ifelse(map$GEN_REF_1 == toupper(consensus), "|", "~")
        match2 <- ifelse(map$GEN_REF_2 == toupper(consensus), "|", "~")
        match1[map$leftRefI]  <- "[" 
        match2[map$rightRefI] <- "]"
        map$text <- cbind(
            map$GEN_REF_1, 
            match1, 
            consensus, 
            match2,
            map$GEN_REF_2
        )
    
    # reference mode (two alignemnts)
    } else if(mode == "Reference Molecule"){
        map$text <- map$text[, which(map$refMolAlns)]
        match1 <- ifelse(map$GEN_REF_1 == toupper(map$text[, 1]), "|", "~")
        match2 <- ifelse(map$GEN_REF_2 == toupper(map$text[, 2]), "|", "~")  
        match1[map$leftRefI]  <- "[" 
        match2[map$rightRefI] <- "]"
        map$text <- cbind(
            map$GEN_REF_1, 
            match1, 
            map$text, 
            match2,
            map$GEN_REF_2
        )  

    # all molecules, all alignments
    } else if(mode == "All Molecules"){
        map$text <- cbind(
            map$GEN_REF_1, 
            map$text, 
            map$GEN_REF_2
        )  
    } else {
        return(NULL)
    }

    # trim the excess
    map$text[is.na(map$text)] <- "~"   
    map$text <- map$text[map$usedPos, ]

    # parse into readable rows
    nChar  <- nrow(map$text)    
    nLines <- ncol(map$text)
    nChunks <- ceiling(nChar / charPerLine)
    nCharLastLine <- nChar %% charPerLine
    if(nCharLastLine == 0) nCharLastLine <- charPerLine
    x <- unlist(sapply(seq_len(nChunks), function(chunk){        
        lineStart <- 1 + (chunk - 1) * charPerLine
        x <- unlist(sapply(seq_len(nLines), function(line){
            bases <- map$text[lineStart:min(nChar, lineStart + charPerLine - 1), line]
            if(any(bases != "~")) {
                bases <- paste(bases, collapse = "")
                bases <- gsub('A', '<span class="base_A">A</span>', bases)
                bases <- gsub('C', '<span class="base_C">C</span>', bases)
                bases <- gsub('G', '<span class="base_G">G</span>', bases)
                bases <- gsub('T', '<span class="base_T">T</span>', bases)
             } else NULL
        }))
        i <- length(x)
        if(i > 2) {
            x[1] <- paste0('<span class="referenceGenome">', x[1], '</span>')
            sapply(2:(i - 2), function(j) x[j] <<- paste0('<span class="alignment">', x[j], '</span>'))
            x[i] <- paste0('<span class="referenceGenome">',  x[i], '</span>')
            paste(x, collapse = "<br>")
        } else NULL
    }))
    x <- if(length(x) > 0) paste(x, collapse = "<br><br>") else NULL
    x <- gsub('~', '&nbsp;', x)

    # annotate the junction positions
    title1 <- paste(map$sv$CHROM_1, map$sv$POS_1, sep = ":")
    title2 <- paste(map$sv$CHROM_2, map$sv$POS_2, sep = ":")
    x <- gsub('[', paste0('<span style="cursor: pointer;" title="', title1, '">*</span>'), x, fixed = TRUE)
    x <- gsub(']', paste0('<span style="cursor: pointer;" title="', title2, '">*</span>'), x, fixed = TRUE)

    return(x)

    # # # app/prepend the outermost chromosomal positions
    # # firstInPair <- 64
    # # svx <- if(bitwAnd(map$sv$FLAG1, firstInPair)) {
    # #     map$sv[,.(RNAME1,PROX_OUT_POS1,RNAME2,PROX_OUT_POS2)]
    # # } else {
    # #     map$sv[,.(RNAME2,PROX_OUT_POS2,RNAME1,PROX_OUT_POS1)]
    # # }
    # # x <- c(
    # #     "<br>",
    # #     paste(svx[[1]], format(svx[[2]], big.mark=","), sep=":"),
    # #     '&darr;',
    # #     x,
    # #     paste(paste(rep("&nbsp;", nCharLastLine-1),collapse=""),'&uarr;',sep=""),
    # #     paste(svx[[3]], format(svx[[4]], big.mark=","), sep=":")
    # # )
    
    # return the hopefully readable sequence block
    paste(x, collapse = "<br>")
}
