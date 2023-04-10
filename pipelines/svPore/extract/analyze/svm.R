# implement support vector machines (SVMs) for splitting chimeric molecules
# train the models using adapters from the outsides of simple alignments 
# logic and workflow:
#   explore the region immediately adjacent to the read outer clips as candidate adapters
#   pad them to allow for sequencing errors and shared bases between adapter and genome
#   use Smith-Waterman local alignment to find best match of clips to expected adapters
#   as training negative controls, use sequences well internal to the clip point, expect to NOT match adapters
#   train an SVM for each end separately, where 3' end clips are generally shorter with fewer bases matched to adapter
#   require that simple alignments have a minimal clip for inclusion in training set
#   provide the SVM with data that summarizes, for paired candidates and controls:
#       the quality of the best match of the sequence to the adapter, e.g., an alignment score
#       the location/position of the match within the sequence
#       together, these allow partial matches but require that they be in a reasonable postion
# use the SVMS to predict adapter sequences at clips, where
#   5'Adapter/----------/3'Adapter...SV_junction...5'Adapter/----------/3'Adapter
#   the alignment query-proximal to the candidate SV junction might clip at a 3' adapter
#   the alignment query-distal   to the candidate SV junction might clip at a 5' adapter
# reject a junction if either adapter matches by removing it from the nodes list

# ONT adapter information
ADAPTER_CORE <- "ACTTCGTTCAGTTACGTATTGCT" # duplex portion of the adapter; last T matches the one-base A-tail
ADAPTER_CORE_RC <- rc(ADAPTER_CORE)      # ADAPTER_CORE is fused to 5' genomic ends, ADAPTER_CORE_RC is fused to 3' ends
coreLen <- nchar(ADAPTER_CORE)
adapterPadding <- 10
outsideLen <- coreLen + adapterPadding
maxCheckLen <- coreLen + 2 * adapterPadding

# SVM shared support functions
getCandidate5 <- function(x){
    if(x$qStart >= outsideLen) list(
        startPos = x$qStart - outsideLen + 1,
        length = maxCheckLen,
        outsideLen = outsideLen
    ) else list(
        startPos = 1,
        length = x$qStart + adapterPadding,
        outsideLen = x$qStart
    )
}
getCandidate3 <- function(x, readLen){
    endPos <- x$qEnd + outsideLen
    overrun <- endPos - readLen
    list(
        startPos = x$qEnd - adapterPadding + 1,
        length = if(overrun <= 0) maxCheckLen else maxCheckLen - overrun
    )
}
runAdapterSW <- function(x, qSeq, ref){
    qry <- substr(qSeq, x$startPos, x$startPos + x$length - 1)
    x$sw <- smith_waterman(qry, ref, fast = FALSE, local = TRUE)
    x
}

# retrieve a sufficient number of 5' and 3' adapter trims for SVM training
extractSvmTrainingSet <- function(nodes, reads){
    message("    extracting a training set of molecule outer clips and random control sequences")
    alns <- reads[molType == edgeTypes$ALIGNMENT]
    d <- do.call(rbind, mclapply(1:nrow(alns), function(i){   
        aln <- alns[i]
        readLen <- nchar(aln$qSeq)
        node <- nodes[aln$qName]

        # 5' side/start of read, expected match to adapter
        # candidate clip sequence
        x5 <- getCandidate5(node)
        x5 <- runAdapterSW(x5, aln$qSeq, ADAPTER_CORE)
        # paired control internal sequence
        x5C <- list(
            startPos = node$qStart + maxCheckLen,
            length = maxCheckLen,
            outsideLen = outsideLen
        )
        x5C <- runAdapterSW(x5C, aln$qSeq, ADAPTER_CORE)

        # 3' side/start of read, expected match to rc(adapter)
        # candidate clip sequence
        x3 <- getCandidate3(node, readLen)
        x3 <- runAdapterSW(x3, aln$qSeq, ADAPTER_CORE_RC)
        # paired control internal sequence
        x3C <- list(
            startPos = node$qEnd - 2 * maxCheckLen,
            length = maxCheckLen
        )
        x3C <- runAdapterSW(x3C, aln$qSeq, ADAPTER_CORE_RC)

        # tabulate
        data.table(
            # 5' side/start of read, expected match to adapter
            # candidate clip sequence
            clip5       = node$qStart, 
            score5      = x5$sw$bestScore, 
            nBases5     = length(x5$sw$qryOnRef),
            start5      = if(x5$sw$bestScore > 0)  x5$sw$qryStart  - x5$outsideLen  - 1 else as.integer(NA),
            end5        = if(x5$sw$bestScore > 0)  x5$sw$qryEnd    - x5$outsideLen  - 1 else as.integer(NA),
            # paired control internal sequence
            score5C     = x5C$sw$bestScore, 
            nBases5C    = length(x5C$sw$qryOnRef),
            start5C     = if(x5C$sw$bestScore > 0) x5C$sw$qryStart - x5C$outsideLen - 1 else as.integer(NA),
            end5C       = if(x5C$sw$bestScore > 0) x5C$sw$qryEnd   - x5C$outsideLen - 1 else as.integer(NA),
            # 3' side/start of read, expected match to rc(adapter)
            # candidate clip sequence
            clip3       = readLen - node$qEnd, 
            score3      = x3$sw$bestScore, 
            nBases3     = length(x3$sw$qryOnRef),
            start3      = if(x3$sw$bestScore > 0)  x3$sw$qryStart  - adapterPadding else as.integer(NA),
            end3        = if(x3$sw$bestScore > 0)  x3$sw$qryEnd    - adapterPadding else as.integer(NA),
            # paired control internal sequence
            score3C     = x3C$sw$bestScore, 
            nBases3C    = length(x3C$sw$qryOnRef),
            start3C     = if(x3C$sw$bestScore > 0) x3C$sw$qryStart - adapterPadding else as.integer(NA),
            end3C       = if(x3C$sw$bestScore > 0) x3C$sw$qryEnd   - adapterPadding else as.integer(NA)
        )   
    }, mc.cores = env$N_CPU))

    # refactor to create the svm table
    rbind(d[, .( # the true clipped ends of reference alignments, for fitting adapter matches
        candidate = TRUE, # controls never have true adapters, candidates usually do, but not always
        trainable5 = clip5 %between% c(15, 50), # whether this clip appears to have an adapter, used to establish the training set
        score5 = score5,
        nBases5 = nBases5,
        start5 = start5,
        end5 = end5,  
        trainable3 = clip3 %between% c(4, 40), # values for the reference clip lengths we train on are set empirically based on plots
        score3 = score3,
        nBases3 = nBases3,
        start3 = start3,
        end3 = end3
    )], d[, .(      # the random paired control sequences for fitting non-matches
        candidate = FALSE,
        trainable5 = TRUE, # all controls are all taken as lacking true adapters and used during training
        score5 = score5C,
        nBases5 = nBases5C,
        start5 = start5C,
        end5 = end5C,
        trainable3 = TRUE,
        score3 = score3C,
        nBases3 = nBases3C,
        start3 = start3C,
        end3 = end3C
    )])
}

# train the svms on all reference and control alignments that have a consistent clip length
# do each side independently
trainAdapterClassifiers <- function(d){
    message("    training the SVM models")
    d <- list(
        "5" = d[, .(candidate = factor(candidate), trainable = trainable5, score = score5, nBases = nBases5, start = start5, end = end5)],
        "3" = d[, .(candidate = factor(candidate), trainable = trainable3, score = score3, nBases = nBases3, start = start3, end = end3)]
    )
    x <- mclapply(d, function(dd){
        svm(candidate ~ score + nBases + start + end, data = dd[trainable == TRUE], scale = TRUE)
    }, mc.cores = 2)
    names(x) <- names(d)
    x
}

# query all SV junctions for possible matches to both adapter sides
extractJunctionSvmParameters <- function(nodes, reads){
    message("    extracting model parameters for all SV junctions")
    I <- nodes[, which(edgeType != edgeTypes$ALIGNMENT)]
    d <- do.call(rbind, mclapply(I, function(i){
        aln3 <- nodes[i - 1] # upstream   of the junction, check for a 3' adapter
        node <- nodes[i]
        aln5 <- nodes[i + 1] # downstream of the junction, check for a 5' adapter
        read <- reads[qName == node[, qName], qSeq] 
        x5 <- getCandidate5(aln5)
        x5 <- runAdapterSW(x5, read, ADAPTER_CORE)  
        x3 <- getCandidate3(aln3, nchar(read))
        x3 <- runAdapterSW(x3, read, ADAPTER_CORE_RC) 
        data.table(
            # unique identifier for this junction (could be more than one per qName)
            qName = node[, qName],
            edge  = node[, edge],
            # downstream of the junction, expected match to adapter
            score5      = x5$sw$bestScore, 
            nBases5     = length(x5$sw$qryOnRef),
            start5      = if(x5$sw$bestScore > 0)  x5$sw$qryStart  - x5$outsideLen  - 1 else as.integer(NA),
            end5        = if(x5$sw$bestScore > 0)  x5$sw$qryEnd    - x5$outsideLen  - 1 else as.integer(NA),
            # upstream of the junction, expected match to rc(adapter)
            score3      = x3$sw$bestScore, 
            nBases3     = length(x3$sw$qryOnRef),
            start3      = if(x3$sw$bestScore > 0)  x3$sw$qryStart  - adapterPadding else as.integer(NA),
            end3        = if(x3$sw$bestScore > 0)  x3$sw$qryEnd    - adapterPadding else as.integer(NA)
        )  
    }, mc.cores = env$N_CPU)) #
}

# use the SVMs to determine whether SV adapter alignments are sufficient evidence to call as an adapter
checkJunctionsForAdapters <- function(svms, d){
    message("    running adapter predictions on SV junctions")
    data.table(
        qName  = d$qName,
        edge   = d$edge,
        score3 = d$score3,
        score5 = d$score5,
        start3 = d$start3,
        end5   = d$end5,
        hasAdapter3 = predict(svms[["3"]], d[, .(score = score3, nBases = nBases3, start = start3, end = end3)]) == TRUE,        
        hasAdapter5 = predict(svms[["5"]], d[, .(score = score5, nBases = nBases5, start = start5, end = end5)]) == TRUE
    )
}

# create two plots for visualizing the results of adapter splitting
plotAdapterMatches <- function(trainingSet, nodes, endN, xcol){
    trainCol    <- paste0("trainable",  endN)
    xcol        <- paste0(xcol,         endN)
    ycol        <- paste0("score",      endN)
    predictCol  <- paste0("hasAdapter", endN)
    trainingSet <- trainingSet[trainingSet[[trainCol]] == TRUE]

    nodes <- nodes[sample(.N, min(500, .N))]

    pngFile <- paste(env$PLOT_PREFIX, "adapters", endN, "png", sep = ".")
    png(pngFile, width = 3, height = 3, units = "in", pointsize = 7, res = 600, type = "cairo")
    plot(
        jitter(trainingSet[[xcol]], amount = 0.5), 
        jitter(trainingSet[[ycol]], amount = 0.5), 
        main = "Training set plus 500 SV junctions",
        xlab = xcol,
        ylab = ycol,
        pch = 19, cex = 0.25, 
        col = "grey"
    )
    points(
        jitter(nodes[[xcol]], amount = 0.5), 
        jitter(nodes[[ycol]], amount = 0.5), 
        pch = 19, cex = 0.25, 
        col = ifelse(nodes[[predictCol]] == TRUE, "red3", "blue")
    )
    dev.off()
}
