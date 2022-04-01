#-------------------------------------------------------------------------------------
# collect outer clips from all samples as additional SV evidence
#-------------------------------------------------------------------------------------
clipFilterScript <- file.path(env$ACTION_DIR, 'find', 'filter_clips.pl')
getOuterClipEvidence <- function(){
    lapply(env$SAMPLES, function(sample){ # work one sample at a time (slow steps parallelized below)
        message(paste(" ", sample))

        # set the input and script files
        compilePrefix <- if(env$FIND_MODE == "compare"){
            prefix <- paste(sample, env$GENOME, 'compile', sep = ".")
            file.path(env$TASK_DIR, sample, prefix)
        } else {
            env$COMPILE_PREFIX
        }
        clipsFile <- paste(compilePrefix, 'outer_clips', 'txt', sep = ".") # TODO: change this to non-indexed gz file

        # load data from the filtering pipe as assembled node molecules
        pipe <- paste("cat", clipsFile, "|", "perl", clipFilterScript, sample, refNodes$file)
        clipMols <- fread(
            cmd = pipe,
            sep = "\t",        
            header = FALSE,
            stringsAsFactors = FALSE,
            col.names = names(compile$working2),
            colClasses = unname(unlist(compile$working2))
        )
        if(nrow(clipMols) == 0) return(NULL)

        # process the clips molecules for duplicates and inversion flipping
        clipMols <- do.call(rbind, mclapply(clipMols[, unique(svIndex)], function(svIdx){ 
            flipInversionClips(svIdx, purgeDuplicateMolecules_(clipMols[svIndex == svIdx]) )
        }, mc.cores = env$N_CPU))
        clipMols
    })
}
#-------------------------------------------------------------------------------------
# flip inversion outer clips as needed to also put them on the canonical strand
#-------------------------------------------------------------------------------------
flipInversionClips <- function(svIdx, clipMols){
    if(!flipGuidance[svIdx, 'isInversion']) return(clipMols)
    flipNode <- unlist(flipGuidance[svIdx, 'flipNode'])
    clipNodes <- clipMols[, ifelse(NODE_1 == "*", 2L, 1L)]
    if(flipNode == 1){
        clipMols[clipNodes == flipNode, ':='(
            SEQ_1      = rc(SEQ_1),
            CLIP_SEQ_1 = rc(CLIP_SEQ_1),
            CIGAR_1    = rc_cigar(CIGAR_1)
        )]
    } else {
        clipMols[clipNodes == flipNode, ':='(
            SEQ_2      = rc(SEQ_2),
            CLIP_SEQ_2 = rc(CLIP_SEQ_2),
            CIGAR_2    = rc_cigar(CIGAR_2)
        )]
    }
    clipMols
}
#-------------------------------------------------------------------------------------
# purge clip molecules
#   from unsequenced, gap-only junctions; they cannot be validated as meaningful
#   that don't match the assembled junction on the unaligned side
#-------------------------------------------------------------------------------------
# assess whether an unaligned sequences, e.g., a clip, matches expectations of aligned reads
getFractionMatchingBases <- function(seq1, seq2, fixedLength = NULL, side = NULL){
    if(is.null(fixedLength)){
        seq1Len <- nchar(seq1)
        seq2Len <- nchar(seq2)
        if(seq1Len > seq2Len){
            seq1 <- if(side == "L") substr(seq1, 1, seq2Len)
                    else substr(seq1, (seq1Len - seq2Len) + 1, seq1Len)
        } else if(seq2Len > seq1Len) {
            seq2 <- if(side == "L") substr(seq2, 1, seq1Len)
                    else substr(seq2, (seq2Len - seq1Len) + 1, seq2Len)
        }
        fixedLength <- nchar(seq1)
    }
    sum(mapply('==', strsplit(seq1, ''), strsplit(seq2, ''))) / fixedLength
}
purgeInvalidClips <- function(seq1, seq2, fixedLength = NULL, side = NULL){
    clipIs <- nodes[, .I[NODE_CLASS == nodeClasses$OUTER_CLIP]]
    if(length(clipIs) > 0){
        discardClipIs <- if(is.na(call$JXN_SEQ)){
            clipIs # result therefore leaves only gap molecules
        } else {
            na.omit(sapply(clipIs, function(i){
                nodeN <- nodes[i, NODE_N]
                x <- getFractionMatchingBases(nodes[i, CLIP_SEQ],
                                              refNodes[nodeN, CLIP_SEQ],
                                              side = refNodes[nodeN, side])
                if(x < env$MIN_MERGE_DENSITY) i else NA
            }))
        }       
        if(length(discardClipIs) > 0) {
            nodes <- nodes[!(MOL_ID %in% nodes[discardClipIs, unique(MOL_ID)])] # remove clips as evidence
        }
    }
}
