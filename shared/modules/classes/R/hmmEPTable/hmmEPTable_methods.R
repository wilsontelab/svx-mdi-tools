#----------------------------------------------------------------------
# hmmEPTable class generic methods, called as method(obj)
#----------------------------------------------------------------------

# solve the HMM for all or a subset of observations to find the most likely path
viterbi <- function(x, ...) {
    UseMethod("viterbi", x)
}
viterbi.hmmEPTable <- function(hmm, observations = TRUE, reverse = FALSE){
    ep <- hmm$emissProbs[observations, ]
    if(reverse) ep <- ep[nrow(ep):1, ]
    tp <- hmm$transProbs
    ep[is.na(ep)] <- -Inf # block unusable paths as having zero probability
    ep[apply(ep, 1, max) == -Inf, ] <- log(1) # mask unusable bins, i.e., those with no usable paths

    # 1. initialization (observation t=1)
    T          <- nrow(ep) # length of the sequence of observations
    N          <- ncol(ep) # number of states
    delta      <- log(matrix(0, nrow = T, ncol = N))
    delta[1, ] <- sapply(1:N, function(i) log(1 / N) + ep[1, i])
    phi        <- matrix(NA, nrow = T, ncol = N)

    # 2. recursion;
    # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
    for (t in 2:T){
        pt <- t - 1
        for (j in 1:N){     # j = this hs
            ep_j <- ep[t, j]
            for (i in 1:N){ # i = prev hs
                delta_ <- delta[pt, i] + tp[i, j] + ep_j
                if(delta[t, j] < delta_){
                    delta[t, j] <- delta_
                    phi[pt, j]  <- i
                }
            }
        }
    }
    
    # 3. termination
    prob <- -Inf
    hsi  <- rep(1, T)
    for (j in 1:N){
        if(prob < delta[T, j]){
            prob <- delta[T, j]
            hsi[T] <- j
        }
    }
    
    # 4. reconstruction and return the hidden state indices
    for (t in (T - 1):1) hsi[t] <- phi[t, hsi[t + 1]]
    if(reverse) rev(hsi) else hsi
}

# run HMM on all unique values of a key vector (e.g., by chromosome)
#----------------------------------------------------------------------
keyedViterbi <- function(x, ...) {
    UseMethod("keyedViterbi", x)
}
keyedViterbi.hmmEPTable <- function(hmm, reverse = FALSE){ # keys is a vector with one value per observation
    if(is.null(hmm$keys)) return( viterbi(hmm, reverse = reverse) )
    ends   <- cumsum(hmm$keys$lengths)
    starts <- c(1, ends + 1)
    unlist(sapply(seq_along(ends), function(i){
        if(starts[i] == ends[i]) return(NA)
        viterbi(hmm, observations = starts[i]:ends[i], reverse = reverse)
    }))
}
