# support functions for analyzeSNVs_server.R

#----------------------------------------------------------------------
# REF/HAP/JXN sequence alignments
#----------------------------------------------------------------------
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

#----------------------------------------------------------------------
# tabulate SNV/indel counts and rates
#----------------------------------------------------------------------
isInformativePos <- function(INF) {
    strsplit(INF, "")[[1]] == "1"
}
isSequencedPos <- function(JXN, MATCH){
    strsplit(JXN, "")[[1]]   != "N" &
    strsplit(MATCH, "")[[1]] != "~"
}
isInterrogatedPos <- function(INF, JXN, MATCH) {
    isInformativePos(INF) & isSequencedPos(JXN, MATCH)
}
isSvVariantPos <- function(MATCH, variantBaseSymbols) {
    strsplit(MATCH, "")[[1]] %in% variantBaseSymbols
}
hapVarSymbols <- c("A", "C", "G", "T", "-", "+")
isSourceVariantPos <- function(HAP1, HAP2){
    strsplit(HAP1, "")[[1]] %in% hapVarSymbols | 
    strsplit(HAP2, "")[[1]] %in% hapVarSymbols
}
tabulateSmallVariants <- function(filteredSvs, settings){
    svs <- filteredSvs()
    req(svs)
    req(nrow(svs) > 0)
    options <- settings$Variant_Options()
    req(options)
    if(options$Allow_Reference_Matches$value) {
        matchThreshold <- SVX$matchTypes$MISMATCH
        variantBaseSymbols <- c("X")
    } else {
        matchThreshold <- SVX$matchTypes$REFERENCE
        variantBaseSymbols <- c("X", ".")
    }
    x <- svs[, .(
        jxnType = svx_jxnType_codeToX(edgeType, "name"),
        hasSmallVariant = MATCH_TYPE >= matchThreshold,
        nInterrogatedBases = {
            sum(isInterrogatedPos(INF_1, JXN_1, MATCH_1)) + 
            sum(isInterrogatedPos(INF_2, JXN_2, MATCH_2))
        },
        nSourceVariantBases = {
            sum(isInterrogatedPos(INF_1, JXN_1, MATCH_1) & isSourceVariantPos(HAP1_1, HAP2_1)) + 
            sum(isInterrogatedPos(INF_2, JXN_2, MATCH_2) & isSourceVariantPos(HAP1_2, HAP2_2)) 
        },
        nSvVariantBases = {
            sum(isInterrogatedPos(INF_1, JXN_1, MATCH_1) & isSvVariantPos(MATCH_1, variantBaseSymbols)) + 
            sum(isInterrogatedPos(INF_2, JXN_2, MATCH_2) & isSvVariantPos(MATCH_2, variantBaseSymbols))
        }
    ), by = SV_ID]
    if(options$Aggregate_Counts$value){
        x <- x[, .(
            nSVs = .N,   
            nSVsWithVariant = sum(hasSmallVariant),               
            nInterrogatedBases = sum(nInterrogatedBases),
            nSourceVariantBases = sum(nSourceVariantBases),
            nSvVariantBases = sum(nSvVariantBases)
        ), by = c("jxnType")]
        x <- rbind(x, x[, .(
            jxnType = "All",
            nSVs = sum(nSVs),
            nSVsWithVariant = sum(nSVsWithVariant), 
            nInterrogatedBases = sum(nInterrogatedBases),
            nSourceVariantBases = sum(nSourceVariantBases),
            nSvVariantBases = sum(nSvVariantBases)
        )])
        x[, ":="(
            svFrequency = nSVsWithVariant / nSVs,
            baseFrequency = nSvVariantBases / nInterrogatedBases
        )]
    }
    x
}

#----------------------------------------------------------------------
# frequency of SNVs occurence by distance from junction
#----------------------------------------------------------------------
distanceSimulation <- function(variant, interrogated, nboots){
    agg <- aggregate(interrogated, list(interrogated), length)
    agg$freq <- agg[[2]] / sum(agg[[2]])
    nVariant <- length(variant)
    medianVariant <- median(variant)
    lessThanMedianVariant <- sapply(1:nboots, function(i){
        x <- sample(agg[[1]], nVariant, replace = TRUE, prob = agg$freq)
        median(x) < medianVariant 
    })
    list("P-Value" = max(sum(lessThanMedianVariant), 1) / nboots)
}
distanceMannWhitney <- function(variant, interrogated, nboots){
    list("P-Value" = wilcox.test(variant, interrogated)$p.value)
}
plotSnvsByDistance <- function(filteredSvs, settings, plot){
    svs <- filteredSvs()
    req(svs)
    svs <- svs[HAS_VARIANT == TRUE]
    req(nrow(svs) > 0)
    options <- settings$Variant_Options()
    req(options)
    startSpinner(session, 'plotSnvsByDistance')

    # collect vectors of interrogated and variant base _distances_ from junction over all SVs in the set
    variantBaseSymbols <- if(options$Allow_Reference_Matches$value) c("X") else c("X", ".")
    interrogated <- svs[, .(distance = c(
        which(rev(isInterrogatedPos(INF_1, JXN_1, MATCH_1))),
        which(    isInterrogatedPos(INF_2, JXN_2, MATCH_2) )
    )), by = SV_ID][, distance]
    variant <- svs[, .(distance = c(
        which(rev(isInterrogatedPos(INF_1, JXN_1, MATCH_1) & isSvVariantPos(MATCH_1, variantBaseSymbols))),
        which(    isInterrogatedPos(INF_2, JXN_2, MATCH_2) & isSvVariantPos(MATCH_2, variantBaseSymbols) )
    )), by = SV_ID][, distance]

    # determine whether the distance distributions appear the same or different
    test <- plot$settings$get("Statistics", "Test_Statistic")
    isTest <- test != "none"
    if(isTest){
        iter <- plot$settings$get("Statistics", "N_Iterations")        
        testFn <- switch(test,
            "none" = NA,
            "Simulated Median"   = distanceSimulation,
            "Mann-Whitney"       = distanceMannWhitney, # for consistency with twosamples     
            "DTS"                = twosamples::dts_test,
            "Anderson-Darling"   = twosamples::ad_test,
            "Cramer-Von Mises"   = twosamples::cvm_test,
            "Kolmogorov-Smirnov" = twosamples::ks_test
        )   
        result <- testFn(variant, interrogated, nboots = iter)
        pValue <- result["P-Value"]

    }

    # plot the two ECDFs
    interrogated <- ecdf(interrogated)
    variant      <- ecdf(variant)
    y <- interrogated(knots(interrogated))    
    stopSpinner(session, 'plotSnvsByDistance')
    plot$initializeFrame(
        xlim = c(0, max(which(y <= 0.99))),
        ylim = c(0, 1),
        xlab = "Distance from Junction (bp)",
        ylab = "Cumulative Frequency",
        title = paste(
            plot$settings$get("Plot_Frame", "Title"),
            if(isTest) paste0(
                "(p <= ", 
                format(pValue, digits = 2, scientific = pValue < 1 / 1000),
                ")"
            ) else ""   
        )  
    )
    plot(
        interrogated, 
        verticals = TRUE, 
        do.points = FALSE, 
        add = TRUE, 
        col = CONSTANTS$plotlyColors$grey
    )
    plot(
        variant, 
        verticals = TRUE, 
        do.points = FALSE, 
        add = TRUE, 
        col = CONSTANTS$plotlyColors$blue
    )
    plot$addLegend(
        legend = c("Sequenced", "Variant"),
        col = c(
            CONSTANTS$plotlyColors$grey,
            CONSTANTS$plotlyColors$blue
        )
    )
}
