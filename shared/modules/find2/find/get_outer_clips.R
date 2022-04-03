#-------------------------------------------------------------------------------------
# collect outer clips from all samples as additional SV evidence
#-------------------------------------------------------------------------------------
clipFilterScript <- file.path(env$ACTION_DIR, 'find', 'filter_clips.pl')
getOuterClipEvidence <- function(){
    lapply(env$SAMPLES, function(sample){ # work one sample at a time (slow steps are parallelized below)
        message(paste(" ", sample))

        # set the input and script files
        compilePrefix <- if(env$FIND_MODE == "compare"){
            prefix <- paste(sample, env$GENOME, 'compile', sep = ".")
            file.path(env$TASK_DIR, sample, prefix)
        } else {
            env$COMPILE_PREFIX
        }
        clipsFile <- paste(compilePrefix, 'outer_clips', 'gz', sep = ".")

        # load data from the filtering pipe as assembled single-node molecules
        pipe <- paste("zcat", clipsFile, "|", "perl", clipFilterScript, sample, refNodesFile)
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
        setkey(clipMols, svIndex)
        do.call(rbind, mclapply(clipMols[, unique(svIndex)], function(svIdx){ 
            flipInversionClips(svIdx, purgeDuplicateMolecules_(clipMols[svIdx]) )
        }, mc.cores = env$N_CPU))
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
purgeInvalidClips <- function(svIdx){

    # check for something to do
    jxnMols <- jxnMols[svIdx] 
    clipIs <- jxnMols[, .I[NODE_CLASS == nodeClasses$OUTER_CLIP]]    
    if(length(clipIs) == 0) return(jxnMols)

    # if not a split or reconstructed gap, purge all clips, their accuracy cannot be checked
    refMol <- jxnMols[IS_REFERENCE == 1]
    discardClipIs <- if(refMol[, JXN_SEQ == "*"]){
        clipIs # result therefore leaves only gap molecules

    # otherwise, check clip sequences against the reference node's clip
    } else {
        na.omit(sapply(clipIs, function(i){
            nodeN <- jxnMols[i, if(NODE_1 == "*") 2 else 1]
            x <- if(nodeN == 1){
                getFractionMatchingBases(
                    jxnMols[i, CLIP_SEQ_1],
                    refMol[  , CLIP_SEQ_1],
                    side = refMol[, side1]
                )
            } else {
                getFractionMatchingBases(
                    jxnMols[i, CLIP_SEQ_2],
                    refMol[  , CLIP_SEQ_2],
                    side = refMol[, side2]
                )
            }
            if(x < env$MIN_MERGE_DENSITY) i else NA
        }))
    }       

    # remove invalid/untrustworthy clips as evidence
    if(length(discardClipIs) > 0) jxnMols <- jxnMols[-discardClipIs] 
    jxnMols
}
