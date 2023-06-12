#-------------------------------------------------------------------------------------
# implement support vector machines (SVMs) for splitting chimeric reads
# train the models using adapters from the outsides of simple, non-SV alignments 
#-------------------------------------------------------------------------------------
# logic and workflow:
#   explore the region immediately adjacent to the read outer clips as candidate adapters
#   pad them to allow for sequencing errors and shared bases between adapter and genome
#   use Smith-Waterman local alignment to find best match of clips to expected adapters
#   as training negative controls, use sequences well internal to the clip point, expect to NOT match adapters
#   train an SVM for each end separately, where 3' end clips are generally shorter with fewer bases matched to adapter
#   require that simple alignments have a minimal clip length for inclusion in the training set
#   provide the SVM with data that summarizes, for paired candidates and controls:
#       the quality of the best match of the sequence to the adapter, e.g., an alignment score
#       the location/position of the match within the sequence
#       together, these allow partial matches but require that they be in a reasonable postion
#   use the SVMs to predict adapter sequences at clips, where
#       5'Adapter/----------/3'Adapter...SV_junction...5'Adapter/----------/3'Adapter
#       the alignment query-proximal to the candidate SV junction might clip at a 3' adapter
#       the alignment query-distal   to the candidate SV junction might clip at a 5' adapter
#   to call a junction as having an adapter, require that it have an SVM hit plus either:
#       a minimal unaligned insertion size, or
#       an unambiguous SW score
#   reject a junction if either adapter matches at or near the unaligned insertion
#-------------------------------------------------------------------------------------

# # ONT adapter information
# # TODO: implement support for mosaic ends in Tn5-based rapid kit
# ADAPTER_CORE <- "ACTTCGTTCAGTTACGTATTGCT" # duplex portion of the adapter; last T matches the one-base A-tail
# ADAPTER_CORE_RC <- rc(ADAPTER_CORE)      # ADAPTER_CORE is fused to 5' genomic ends, ADAPTER_CORE_RC is fused to 3' ends
# coreLen <- nchar(ADAPTER_CORE)
# adapterPadding <- 10
# outsideLen <- coreLen + adapterPadding
# maxCheckLen <- coreLen + 2 * adapterPadding

# # SVM shared support functions
# getCandidate5 <- function(edge){
#     if(edge$qStart >= outsideLen) list(
#         startPos = edge$qStart - outsideLen + 1,
#         length = maxCheckLen,
#         outsideLen = outsideLen
#     ) else list(
#         startPos = 1,
#         length = edge$qStart + adapterPadding,
#         outsideLen = edge$qStart
#     )
# }
# getCandidate3 <- function(edge, readLen){
#     endPos <- edge$qEnd + outsideLen
#     overrun <- endPos - readLen
#     list(
#         startPos = edge$qEnd - adapterPadding + 1,
#         length = if(overrun <= 0) maxCheckLen else maxCheckLen - overrun
#     )
# }
# runAdapterSW <- function(x, qSeq, ref){
#     qry <- substr(qSeq, x$startPos, x$startPos + x$length - 1)
#     x$sw <- smith_waterman(qry, ref, fast = FALSE, local = TRUE)
#     x
# }

# retrieve a sufficient number of 5' and 3' adapter trims for SVM training
extractSvmTrainingSet <- function(trainingEdges){
# extractSvmTrainingSet <- function(reads, trainingEdges){
#     message("extracting a training set of read outer clips and random control sequences")
#     edges <- trainingEdges[qEnd - qStart > 300 & gapCompressedIdentity > env$MIN_ALIGNMENT_IDENTITY]    
#     reads <- reads[qName %in% edges[, qName]]
#     setkey(edges, qName)
#     d <- do.call(rbind, mclapply(1:nrow(reads), function(i){   
#         read <- reads[i]
#         edge <- edges[read[, qName]]
#         readLen <- nchar(read$qSeq)

#         # 5' side/start of read, expected match to adapter
#         # candidate clip sequence
#         x5 <- getCandidate5(edge)
#         x5 <- runAdapterSW(x5, read$qSeq, ADAPTER_CORE)
#         # paired control internal sequence
#         x5C <- list(
#             startPos = edge$qStart + maxCheckLen,
#             length = maxCheckLen,
#             outsideLen = outsideLen
#         )
#         x5C <- runAdapterSW(x5C, read$qSeq, ADAPTER_CORE)

#         # 3' side/start of read, expected match to rc(adapter)
#         # candidate clip sequence
#         x3 <- getCandidate3(edge, readLen)
#         x3 <- runAdapterSW(x3, read$qSeq, ADAPTER_CORE_RC)
#         # paired control internal sequence
#         x3C <- list(
#             startPos = edge$qEnd - 2 * maxCheckLen,
#             length = maxCheckLen
#         )
#         x3C <- runAdapterSW(x3C, read$qSeq, ADAPTER_CORE_RC)

#         # tabulate
#         x <- data.table(
#             # 5' side/start of read, expected match to adapter
#             # candidate clip sequence
#             clip5       = edge$qStart, 
#             score5      = x5$sw$bestScore, 
#             nBases5     = length(x5$sw$qryOnRef),
#             start5      = if(x5$sw$bestScore > 0)  x5$sw$qryStart  - x5$outsideLen  - 1 else as.integer(NA),
#             end5        = if(x5$sw$bestScore > 0)  x5$sw$qryEnd    - x5$outsideLen  - 1 else as.integer(NA),
#             # paired control internal sequence
#             score5C     = x5C$sw$bestScore, 
#             nBases5C    = length(x5C$sw$qryOnRef),
#             start5C     = if(x5C$sw$bestScore > 0) x5C$sw$qryStart - x5C$outsideLen - 1 else as.integer(NA),
#             end5C       = if(x5C$sw$bestScore > 0) x5C$sw$qryEnd   - x5C$outsideLen - 1 else as.integer(NA),
#             # 3' side/start of read, expected match to rc(adapter)
#             # candidate clip sequence
#             clip3       = readLen - edge$qEnd, 
#             score3      = x3$sw$bestScore, 
#             nBases3     = length(x3$sw$qryOnRef),
#             start3      = if(x3$sw$bestScore > 0)  x3$sw$qryStart  - adapterPadding else as.integer(NA),
#             end3        = if(x3$sw$bestScore > 0)  x3$sw$qryEnd    - adapterPadding else as.integer(NA),
#             # paired control internal sequence
#             score3C     = x3C$sw$bestScore, 
#             nBases3C    = length(x3C$sw$qryOnRef),
#             start3C     = if(x3C$sw$bestScore > 0) x3C$sw$qryStart - adapterPadding else as.integer(NA),
#             end3C       = if(x3C$sw$bestScore > 0) x3C$sw$qryEnd   - adapterPadding else as.integer(NA)
#         )  
#     }, mc.cores = env$N_CPU))

    # refactor to create the svm table
    rbind(trainingEdges[, .( # the true clipped ends of reference alignments, for fitting adapter matches
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
    )], trainingEdges[, .(      # the random paired control sequences for fitting non-matches
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
    message("training the SVM models")
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

# # query all SV junctions for possible matches to both adapter sides
# extractJunctionSvmParameters <- function(edges, reads){
#     message("extracting model parameters for all SV junctions")
#     I <- edges[, which(edgeType != edgeTypes$ALIGNMENT)]
#     setkey(reads, qName)
#     d <- do.call(rbind, mclapply(I, function(i){
#         aln3 <- edges[i - 1] # upstream   of the junction, check for a 3' adapter
#         edge <- edges[i]
#         aln5 <- edges[i + 1] # downstream of the junction, check for a 5' adapter
#         qSeq <- reads[edge[, qName], qSeq]
#         x5 <- getCandidate5(aln5)
#         x5 <- runAdapterSW(x5, qSeq, ADAPTER_CORE)  
#         x3 <- getCandidate3(aln3, nchar(qSeq))
#         x3 <- runAdapterSW(x3, qSeq, ADAPTER_CORE_RC) 
#         data.table(
#             # unique identifier for this junction (could be more than one per qName)
#             qName      = edge[, qName],
#             edgeN      = edge[, edgeN],
#             insertSize = edge[, insertSize],
#             # downstream of the junction, expected match to adapter
#             score5      = x5$sw$bestScore, 
#             nBases5     = length(x5$sw$qryOnRef),
#             start5      = if(x5$sw$bestScore > 0)  x5$sw$qryStart  - x5$outsideLen  - 1 else as.integer(NA),
#             end5        = if(x5$sw$bestScore > 0)  x5$sw$qryEnd    - x5$outsideLen  - 1 else as.integer(NA),
#             # upstream of the junction, expected match to rc(adapter)
#             score3      = x3$sw$bestScore, 
#             nBases3     = length(x3$sw$qryOnRef),
#             start3      = if(x3$sw$bestScore > 0)  x3$sw$qryStart  - adapterPadding else as.integer(NA),
#             end3        = if(x3$sw$bestScore > 0)  x3$sw$qryEnd    - adapterPadding else as.integer(NA)
#         )
#     }, mc.cores = env$N_CPU))
# }

# use the SVMs to determine whether SV adapter alignments are sufficient evidence to call as an adapter
checkJunctionsForAdapters <- function(svms, d){
    message("running adapter predictions on SV junctions")
    predictAdapter <- function(endN, insertSize, dd){
        scoreThreshold <- if(endN == "3") 7 else 8 # determined empirically based on ligation kit adapters
        predict(svms[[endN]], dd) == TRUE & (insertSize >= 5 | dd$score > scoreThreshold)
    }
    data.table(
        qName  = d$qName,
        edgeN  = d$edgeN,
        # score3 = d$score3,
        # score5 = d$score5,
        # start3 = d$start3,
        # end5   = d$end5,
        hasAdapter3 = predictAdapter("3", d$insertSize, d[, .(score = score3, nBases = nBases3, start = start3, end = end3)]),        
        hasAdapter5 = predictAdapter("5", d$insertSize, d[, .(score = score5, nBases = nBases5, start = start5, end = end5)])
    )
}

# update edges with adapter splitting
updateEdgesForAdapters <- function(edges, adapterCheck){
    message("updating edges")
    dropCols <- names(adapterCheck)[!(names(adapterCheck) %in% c("qName","edgeN"))]
    edges <- edges[, .SD, .SDcols = names(edges)[!(names(edges) %in% dropCols)]]
    edges <- merge(edges, adapterCheck, by = c("qName","edgeN"), all.x = TRUE)
    edges[, hasAdapter := hasAdapter3 | hasAdapter5]
    setkey(edges, qName, blockN, edgeN)
}
