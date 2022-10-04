#----------------------------------------------------------------------
# object class for constructing and solving a Hidden Markov Model where:
#   emission probs are given as matrix [observation,state], with values as log(ep)
#   transition probs are derived from a single non-log value, with preference to staying in state
#   start probs are random, i.e., equal for all hidden states
#----------------------------------------------------------------------
# providing an equal length vector of keys supports calls to keyedViterbi (e.g., by chrom)
#----------------------------------------------------------------------
# providing an equal length vector of predetermined state Is restricts the HMM output at those observations
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN CLASS
#----------------------------------------------------------------------
new_hmmEPTable <- function(emissProbs, transProb = 1e-6, keys = NULL, forceStateIs = NULL){ 
    class <- 'hmmEPTable' # for reportProgress tracing
    
#----------------------------------------------------------------------
# parse input values
#----------------------------------------------------------------------
T <- nrow(emissProbs)
N <- ncol(emissProbs)
if(!is.null(forceStateIs)){ # override emissProbs for predetermined states
    forcedProbs <- lapply(1:max(forceStateIs, na.rm = TRUE), function(i){
        ep <- log( rep(0, N) )
        ep[i] <- log(1)
        ep
    })
    emissProbs <- t(sapply(1:T, function(i){
        j <- forceStateIs[i]
        if(is.na(j)) emissProbs[i, ] else forcedProbs[[j]]
    }))  
}
transProbs <- matrix(transProb, nrow = N, ncol = N)
diag(transProbs) <- 1 - transProb * (N - 1)
transProbs <- log(transProbs)
keys <- if(is.null(keys)) NULL else rle(keys)

#----------------------------------------------------------------------
# set the return value
#----------------------------------------------------------------------
structure(
    list(
        emissProbs = emissProbs,
        transProbs = transProbs,
        keys = keys
    ),
    class = class
)

#----------------------------------------------------------------------
# END CLASS
#----------------------------------------------------------------------
}
#----------------------------------------------------------------------
