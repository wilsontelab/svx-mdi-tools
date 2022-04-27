#----------------------------------------------------------------------
# functions to load SVs and associates evidence molecules
#----------------------------------------------------------------------

# ----------------------------------------------------------------------
# map of all base positions at and around the selected SV junction
# supports consensus making, evident plotting, etc.
# ----------------------------------------------------------------------
wsh <- 0.7 # for washed-out colors
zoomBaseColors <- list( # generally follow IGV base color conventions
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
clipBases <- c('a', 'c', 'g', 't', 'n')
# zoomBaseR <- sapply(zoomBaseColors, function(v) v[1])
# zoomBaseG <- sapply(zoomBaseColors, function(v) v[2])
# zoomBaseB <- sapply(zoomBaseColors, function(v) v[3])

# convert an SV and its set of evidence molecules to a base map
getJunctionMap <- function(x){
    req(x)
    req(x$sv$JXN_BASES != "*")
    startSpinner(session, 'getJunctionMap') 

    ##########
    x$mols <- x$mols[IS_REFERENCE == 1]

    # initialize map and dimensions
    refWidth <- nchar(x$sv$GEN_REF_1)
    faidxPadding <- (refWidth - 1) / 2
    leftRefI <- faidxPadding + 1
    microhomologyLength <- x$sv[, MICROHOM_LEN]
    mapWidth <- refWidth + 1 - microhomologyLength
    nAln <- x$mols[, .N + sum(!isOuterClip)]
    vPadding <- 3
    yOffset <- 1
    baseMap <- matrix(NA, nrow = mapWidth, ncol = nAln) # for consensus making
    imgMap  <- array(NA, dim = c(mapWidth, nAln + vPadding * 2, 3)) # for graphical display

    # parse the genome reference lines
    GEN_REF_1 <- c(    strsplit(x$sv$GEN_REF_1, "")[[1]],  rep(" ", mapWidth))
    GEN_REF_2 <- c(rev(strsplit(x$sv$GEN_REF_2, "")[[1]]), rep(" ", mapWidth))
    GEN_REF_1 <- GEN_REF_1[1:mapWidth]
    GEN_REF_2 <- GEN_REF_2[1:mapWidth]
    GEN_REF_2 <- rev(GEN_REF_2)

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
    processNode <- function(nodeN, jxnSide, jxnPos, CIGAR, SEQ, pos){
        if(CIGAR == "*") return(NULL) # the missing node of an outer clip
        aln <- cigarToRefAln(CIGAR, SEQ)            
        #      if(clipMode == "grey") aln$seq[aln$seq %in% clipBases] <- "x"
        # else if(clipMode == "none") aln$seq[aln$seq %in% clipBases] <- " "

        # aln$seq[aln$seq %in% clipBases] <- " "

        nodePosDelta <- if(jxnSide == "L") jxnPos - pos else pos - jxnPos
        mapPos <- if(nodeN == 1){
            leftRefI - nodePosDelta - aln$lengthAlignedOut - aln$leftClip + 1
        } else {
            leftRefI + 2 - microhomologyLength + nodePosDelta - aln$leftClip
        }
        j <- mapPos:(mapPos + aln$lengthOut - 1)
        col <- t(sapply(aln$seq, function(base) zoomBaseColors[[base]]))
        imgMap[ j, yOffset + vPadding, ] <<- col
        baseMap[j, yOffset] <<- aln$seq
        yOffset <<- yOffset + 1 
    }
    x$mols[
        order( CLIP_LEN_1, POS_1), 
        mapply(processNode, 1, x$sv$SIDE_1, x$sv$POS_1, CIGAR_1, SEQ_1, pos1)
    ]
    x$mols[
        order(-CLIP_LEN_2, POS_2 + nchar(SEQ_2)), 
        mapply(processNode, 2, x$sv$SIDE_2, x$sv$POS_2, CIGAR_2, SEQ_2, pos2)
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
        usedPos = usedPos,
        GEN_REF_1 = GEN_REF_1,
        GEN_REF_2 = GEN_REF_2
    )
}

# create an alignment betweeen junction molecules, their consensus and the genome reference
getJunctionAlignment <- function(map){
    req(map)    

    # consensus <- apply(map$text, 1, function(x){
    #     x <- x[!(is.na(x) | x %in% clipBases)]
    #     if(length(x) == 0) return(NA)
    #     agg <- aggregate(x, list(x), length)
    #     agg[which.max(agg[[2]]), 1]
    # })
    # match1 <- ifelse(map$GEN_REF_1 == consensus, "|", "~")
    # match2 <- ifelse(map$GEN_REF_2 == consensus, "|", "~")

    # map$text <- cbind(
    #     map$GEN_REF_1, 
    #     match1,
    #     consensus, 
    #     match2,        
    #     map$GEN_REF_2, 
    #     rep(" ", length(map$GEN_REF_1))
    # )[map$usedPos, ]

    map$text[is.na(map$text)] <- "~"   

    map$text <- map$text[map$usedPos, ]

    # dstr(map) 

    # parse into readable rows
    nChar  <- nrow(map$text)    
    nLines <- ncol(map$text)
    charPerLine <- 100 # as.numeric(input$charPerLine)
    nChunks <- ceiling(nChar / charPerLine)
    nCharLastLine <- nChar %% charPerLine
    if(nCharLastLine == 0) nCharLastLine <- charPerLine

    # output <- character()
    # blueLines <- list( # junction lines to denote with highlight coloring
    #     Genome = 3:4,
    #     'Genome+' = 3:4,
    #     Reads = 1  
    # )

    # dprint(map$sv[, .(
    #     MERGE_LEN, MICROHOM_LEN, JXN_BASES
    # )])
    # dprint(map$mols[, .(
    #     CIGAR_1, SEQ_1, CLIP_LEN_1, 
    #     CIGAR_2, SEQ_2, CLIP_LEN_2)]
    # )

    return( paste(sapply(seq_len(nChunks), function(chunk){
        # if(chunk > 1) output <- c(output, "<br>")        
        lineStart <- 1 + (chunk - 1) * charPerLine
        paste(sapply(seq_len(nLines), function(line){
            bases <- map$text[lineStart:min(nChar, lineStart + charPerLine - 1), line]
            gsub('~', '&nbsp;', paste(bases, collapse = ""))
        }), collapse = "<br>")
    }), collapse = "<br><br>") )

    # for(chunk in seq_along(nChunks)){
    #     if(chunk > 1) output <- c(output, "<br>")        
    #     lineStart <- 1 + (chunk - 1) * charPerLine
    #     for(line in seq_along(nLines)){
    #         bases <- map$bases[lineStart:(lineStart + charPerLine - 1), line]
    #         str  <- substr(map$bases[j], start, start + charPerLine - 1)

    #         str  <- gsub(' ', '&nbsp;', str) 
    #         str_ <- gsub('~', '', str) 
    #         if(nchar(str_) > 0){
    #             col <- if(j %in% blueLines[[input$jxnDisplayType]]) 'blue' else 'black'
    #             output <- c(
    #                 output,
    #                 paste(paste('<span style="color:', col, ';">', sep = ""),
    #                       str,
    #                       '</span>', sep = "")
    #             )                
    #         }            
    #     }
    # }
    # output <- gsub('~', '&nbsp;', output)
    
    # # # app/prepend the outermost chromosomal positions
    # # firstInPair <- 64
    # # svx <- if(bitwAnd(map$sv$FLAG1, firstInPair)) {
    # #     map$sv[,.(RNAME1,PROX_OUT_POS1,RNAME2,PROX_OUT_POS2)]
    # # } else {
    # #     map$sv[,.(RNAME2,PROX_OUT_POS2,RNAME1,PROX_OUT_POS1)]
    # # }
    # # output <- c(
    # #     "<br>",
    # #     paste(svx[[1]], format(svx[[2]], big.mark=","), sep=":"),
    # #     '&darr;',
    # #     output,
    # #     paste(paste(rep("&nbsp;", nCharLastLine-1),collapse=""),'&uarr;',sep=""),
    # #     paste(svx[[3]], format(svx[[4]], big.mark=","), sep=":")
    # # )
    
    # return the hopefully readable sequence block
    paste(output, collapse = "<br>")
}