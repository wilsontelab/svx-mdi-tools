# support functions for analyzeSNVs_server.R
alnCols1 <-     c("INF_1", "REF_1", "HAP1_1", "HAP2_1", "MATCH_1", "JXN_1")
alnCols2 <- rev(c("INF_2", "REF_2", "HAP1_2", "HAP2_2", "MATCH_2", "JXN_2"))
refCols  <- c("REF", "HAP1", "HAP2")
parseSide <- function(sv, sideI, charPerLine){

    # parse the columns for this junction side
    alnCols <- if(sideI == 1) alnCols1 else alnCols2
    alnColI <- as.list(seq_along(alnCols))
    names(alnColI) <- gsub(paste0("_", sideI), "", alnCols)

    # split the display lines to vectors for parsing
    # complement the side of the junction inside an inverted segment
    # (the genotype action provides all sequences as reversed but not complemented from reference)
    text <- sapply(alnCols, function(col) strsplit(sv[[col]], "")[[1]])
    if(sv$SIDE_1 == sv$SIDE_2){
        flipSide <- if(sv$SIDE_1 == "R") 1 else 2
        if(sideI == flipSide) for(col in c(refCols, "JXN")){
            text[, alnColI[[col]]] <- rc(text[, alnColI[[col]]])
        }
    }
    
    # parse the display
    nChar <- nrow(text)   
    nLines <- ncol(text)
    nChunks <- ceiling(nChar / charPerLine)
    nCharLastLine <- nChar %% charPerLine
    if(nCharLastLine == 0) nCharLastLine <- charPerLine
    x <- unlist(sapply(seq_len(nChunks), function(chunk){        
        lineStart <- 1 + (chunk - 1) * charPerLine
        x <- unlist(sapply(seq_len(nLines), function(line){
            bases <- text[lineStart:min(nChar, lineStart + charPerLine - 1), line]
            if(any(bases != "~")) {
                bases <- paste(bases, collapse = "")
                bases <- gsub('A', '<span class="base_A">A</span>', bases)
                bases <- gsub('C', '<span class="base_C">C</span>', bases)
                bases <- gsub('G', '<span class="base_G">G</span>', bases)
                bases <- gsub('T', '<span class="base_T">T</span>', bases)
                bases <- gsub('X', '<span class="mismatch">X</span>', bases)
             } else NULL
        }))
        for(col in refCols) {
            x[alnColI[[col]]] <- paste0('<span class="referenceGenome">', x[alnColI[[col]]], '</span>')
        }
        x[alnColI$JXN] <- paste0('<span class="alignment">', x[alnColI$JXN], '</span>')
        x[alnColI$MATCH] <- paste0('<span class="matches">', x[alnColI$MATCH], '</span>')
        x[alnColI$INF] <- gsub("1", "=", x[alnColI$INF])
        x[alnColI$INF] <- gsub("0", " ", x[alnColI$INF])
        paste(x, collapse = "<br>")
    }))
    x <- if(length(x) > 0) paste(x, collapse = "<br><br>") else NULL
    x <- gsub('~', '&nbsp;', x)
    x <- gsub('-', '&#8209;', x)
    x
}

# tabulate SNV/indel counts and rates
informativeSeqCols <- c("INF_1",   "INF_2")
informativeSeqCols <- c("MATCH_1", "MATCH_1")

alnCols1 <-     c("", "REF_1", "HAP1_1", "HAP2_1", "MATCH_1", "JXN_1")
alnCols2 <- rev(c("", "REF_2", "HAP1_2", "HAP2_2", "MATCH_1", "JXN_2"))

tabulateSmallVariants <- function(filteredSvs, settings){
    svs <- filteredSvs()
    req(svs)
    svFilters <- settings$Variant_Filters()
    req(svFilters)
    if(svFilters$Allow_Reference_Matches$value) {
        matchThreshold <- SVX$matchTypes$MISMATCH
        variantBaseSymbols <- c("X")
    } else {
        matchThreshold <- SVX$matchTypes$REFERENCE
        variantBaseSymbols <- c("X", ".")
    }
    setkey(SVX$jxnTypes, code)   



    svs[, .(
        JXN_TYPE = SVX$jxnTypes[JXN_TYPE, name],
        hasSmallVariant = ifelse(MATCH_TYPE >= matchThreshold, TRUE, FALSE),
        nBasesInterrogated = {
            sum(
                strsplit(MATCH_1, "")[[1]] != "~" & 
                strsplit(INF_1, "")[[1]]   == "1"
            ) + 
            sum(
                strsplit(MATCH_2, "")[[1]] != "~" & 
                strsplit(INF_2, "")[[1]]   == "1"
            )
        },
        nVariantBases = {
            sum(
                strsplit(MATCH_1, "")[[1]] %in% variantBaseSymbols & 
                strsplit(INF_1, "")[[1]]   == "1"
            ) + 
            sum(
                strsplit(MATCH_2, "")[[1]] %in% variantBaseSymbols & 
                strsplit(INF_2, "")[[1]]   == "1"
            )
        }
    ), by = SV_ID]
}
