
# smith_waterman performs Smith-Waterman alignment of a query to a reference sequence
# see smith_waterman.pl for more description

# constants for matrix traversal
DIAG_  <- 1
UP_    <- 2
LEFT_  <- 3
QRY_START <- 1
QRY_END   <- 2

# set scoring metrics to user provided values or defaults
matchScore          <- 1
mismatchPenalty     <- -1.5
gapOpenPenalty      <- -2.501 # 0.001 ajustment gives slight preference to not opening a single-base terminal gap
gapExtensionPenalty <- -1
maxShift            <- 3

# create the lookup table for match/mismatch score for all possible IUPAC code combinations
initializePairScores <- function() {
    scores <- list()
    for (key in names(baseMatches))   scores[[key]] <- matchScore      # e.g. A:A, full match
    for (key in names(ryswkmMatches)) scores[[key]] <- matchScore / 2  # e.g. A:R, half match
    for (key in names(nMatches))      scores[[key]] <- 0  # e.g. A:N, uninformative base position, neither promoted nor penalized # nolint
    scores[['--']] <- 0
    bases <- c('A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'N', '-')
    for (base1 in bases){
        for (base2 in bases){
            key <- paste0(base1, base2)
            if(is.null(scores[[key]])) scores[[key]] <- mismatchPenalty
        }
    }
    scores
}
pairedBaseScores <- initializePairScores()

# the Smith-Waterman algorithm itself
smith_waterman <- function(qry, ref, fast = TRUE, forceQryEnd = NULL){ # setting fast to TRUE enforces register shift limitations # nolint

    # collect sequence inputs
    qry <- strsplit(qry, '')[[1]]
    ref <- strsplit(ref, '')[[1]]
    nQ <- length(qry)
    nR <- length(ref)
    isforceQryEnd <- !is.null(forceQryEnd)
    if(isforceQryEnd && forceQryEnd == QRY_START){
        qry <- rev(qry) # temporarily reverse sequence to allow same code for forceQryEnd QRY_START and QRY_END
        ref <- rev(ref) # algorithm is written for QRY_END (i.e. the right side of qry)
    }

    # fill score matrix based on matches and gaps
    x <- array(0, dim = c(nR + 1, nQ + 1, 2))
    bestScore <- 0
    bestPath <- NULL
    paths <- list()
    for(refI in 1:nR){
        refBase <- ref[refI]    
        isTermI <- refI == nR
        minQryI <- if(fast) refI - maxShift else 1
        maxQryI <- if(fast) refI + maxShift else nQ
        for(qryI in minQryI:maxQryI) { # limit the possible query to reference base matches
            if(qryI < 1 || qryI > nQ) next
            diag_score <- x[refI,  qryI,     1] + pairedBaseScores[[ paste0(refBase, qry[qryI]) ]]  
            up_score   <- x[refI,  qryI + 1, 1] + (if(x[refI,  qryI + 1, 2] == UP_)   gapExtensionPenalty else gapOpenPenalty) # nolint
            left_score <- x[refI + 1, qryI,  1] + (if(x[refI + 1, qryI,  2] == LEFT_) gapExtensionPenalty else gapOpenPenalty) # nolint
            if(diag_score >= up_score &&
               diag_score >= left_score){
                score   <- diag_score
                pointer <- DIAG_
            } else if(up_score >= diag_score &&
                      up_score >= left_score){
                score   <- up_score
                pointer <- UP_  
            } else {
                score   <- left_score
                pointer <- LEFT_  
            }       
            if(isforceQryEnd){
                if(score > 0){
                    x[refI + 1, qryI + 1, 1] <- score
                    x[refI + 1, qryI + 1, 2] <- pointer
                }
                if(qryI == nQ){ # ensure that all reported alignments go to end of query
                    if(score > bestScore){
                        bestScore <- score
                        paths <- list(c(refI, qryI))
                    } else if(score == bestScore){
                        paths <- c(paths, c(refI, qryI))
                    }
                }            
            } else { # general untrimmed alignment when requiring end-to-end alignment, e.g. in consensus building
                x[refI + 1, qryI + 1, 1] <- score
                x[refI + 1, qryI + 1, 2] <- pointer
                if((isTermI || qryI == nQ) && score > bestScore){
                    bestScore <- score
                    bestPath  <- c(refI, qryI)
                }                
            }   
        }
    }            

    # trace backwards to deconvolute best matching path(s) and alignment map(s)
    if(isforceQryEnd){ # demand just one best hit in forceQryEnd mode
        if(length(paths) > 1) return(list(score = 0))
        bestPath <- paths[[1]]
    }
    if(is.null(bestPath)) return(list(score = 0))
    maxRefI <- bestPath[1]
    maxQryI <- bestPath[2]
    qryOnRef <- c()
    refI <- maxRefI + 1
    qryI <- maxQryI + 1
    while (1) {
        pointer <- x[refI, qryI, 2]
        if(pointer == 0) break
        if (pointer == DIAG_) { # M operations
            qryOnRef <- c(qry[qryI - 1], qryOnRef)
            refI <- refI - 1
            qryI <- qryI - 1
        } else if (pointer == LEFT_) { # I operation relative to reference, append to NEXT reference base
            qryOnRef[1] <- c(paste0( qry[qryI - 1], if(length(qryOnRef) == 0) "" else qryOnRef[1] ))
            qryI <- qryI - 1
        } else {  # up <- D operation relative to reference, pad with a dummy character
            qryOnRef <- c('-', qryOnRef)
            refI <- refI - 1
        }
    }

    # return best alignment
    if(isforceQryEnd && forceQryEnd == QRY_START){
        qryOnRef <- rev(qryOnRef); # revert back to original sequence orientation when QUERY_START
        tmp <- qryI
        qryI <- nQ - maxQryI + 1
        maxQryI <- nQ - tmp + 1
        tmp <- refI
        refI <- nR - maxRefI + 1
        maxRefI <- nR - tmp + 1
    } 
    return(list(
        qryOnRef = qryOnRef,
        bestScore = bestScore,
        qryStart = qryI,
        qryEnd = maxQryI,
        refStart = refI,
        refEnd = maxRefI       
    ))      
}
