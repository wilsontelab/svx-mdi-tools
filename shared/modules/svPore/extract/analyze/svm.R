#-------------------------------------------------------------------------------------
# implement support vector machines (SVMs) for splitting chimeric reads
# train the models using adapters from the outsides of simple, non-SV alignments 
#-------------------------------------------------------------------------------------
# logic and workflow:
#   explore the region immediately adjacent to read outer clips as candidate adapters
#   pad them to allow for sequencing errors and shared bases between adapter and genome
#   use Smith-Waterman local alignment to find best match of clips to expected adapters (done by extend_edges.pl)
#   as training negative controls, use sequences well internal to the clip point, expected to NOT match adapters
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
# at this point, initial scoring actions are executed by extract/extend_edges.pl
# this R script uses the scores to create and apply the SVMs
#-------------------------------------------------------------------------------------

# refactor svm scores to create the table suitable for training
extractSvmTrainingSet <- function(trainingEdges){
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

# use the SVMs to determine whether SV adapter alignments are sufficient evidence to call a junction as an adapter chimera
checkJunctionsForAdapters <- function(svms, d){
    message("running adapter predictions on SV junctions")
    predictAdapter <- function(endN, insertSize, dd){
        scoreThreshold <- if(endN == "3") 7 else 8 # determined empirically based on ligation kit adapters
        predict(svms[[endN]], dd) == TRUE & (insertSize >= 5 | dd$score > scoreThreshold)
    }
    data.table(
        readI = d$readI,
        edgeN    = d$edgeN,
        hasAdapter5 = predictAdapter("5", d$insertSize, d[, .(score = score5, nBases = nBases5, start = start5, end = end5)]),        
        hasAdapter3 = predictAdapter("3", d$insertSize, d[, .(score = score3, nBases = nBases3, start = start3, end = end3)])        
    )
}

# update edges with adapter splitting flags
# no edge are dropped yet, at this point we only set hasAdapter5 and hasAdapter3
updateEdgesForAdapters <- function(edges, adapterCheck = NULL){
    message("updating edges")
    if(is.null(adapterCheck)){ # SKIP_ADAPTER_CHECK
        edges[edgeType != edgeTypes$ALIGNMENT, ":="(
            hasAdapter5 = FALSE,
            hasAdapter3 = FALSE
        )]
    } else {
        edges <- merge(edges, adapterCheck, by = c("readI","edgeN"), all.x = TRUE)
        edges[, c(edgeAdapterScores) := NULL]
        edges[edgeType != edgeTypes$ALIGNMENT & is.na(hasAdapter5), ":="(
            hasAdapter5 = FALSE, # these are junctions with short or no insertions for which SW+SVM was skipped
            hasAdapter3 = FALSE
        )]
    }
    setkey(edges, readI, blockN, edgeN)
}
