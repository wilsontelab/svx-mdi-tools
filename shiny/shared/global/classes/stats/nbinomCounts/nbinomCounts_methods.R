#----------------------------------------------------------------------
# nbinomCounts class generic methods, called as method(obj)
#----------------------------------------------------------------------

# predict the fit value, analogous to predict.lm, etc.
predict.nbinomCounts <- function(object, newdata=NULL, mu=NULL, width=100, ...){
    if(is.null(newdata)) newdata <- object$data    
    if(is.null(mu)) mu <- object$mu
    #if(method[1] == 'median'){
    #    
    #    
    #    
    #    
    #} else {
        sapply(newdata$x, function(x){ # fit an expected read count for all _requested_ points ...
            is <- order(abs(object$data$x - x))[1:width] # ... based on indices of N _modeled_ points with the closest predictor values # nolint
            new_nbinomCounts( data = object$data[is, ] )$peak        
        })        
    #}
}


# create an HMM object for solving by HHH::viterbi
initHMM.nbinomCounts <- function(nb, # the object created by new_nbinomCounts, which established theta for the data set
                                 data, # an ordered vector of all counts to be modeled
                                 startProbs = NULL, # vector of starting probabilities; defaults to equal for all hidden states # nolint
                                 transProb = 1e-6, # a single, fixed probablity of changing states; should be rare
                                 mu = NULL, # a vector of mu values to be used as the hidden states
                                 modeledCN = NULL, minCN = 0, maxCN = 6 # alternatively to mu, assume a copy number model over this range # nolint
){
    if(is.null(mu)) mu <- nb$mu * minCN:maxCN / modeledCN
    nHiddenStates <- length(mu)
    data <- nb$nonNegInt(data)
    allCounts <- 0:max(data)
    if(is.null(startProbs)) startProbs <- rep(1 / nHiddenStates, nHiddenStates)
    transProbs <- matrix(transProb, nrow = nHiddenStates, ncol = nHiddenStates)
    diag(transProbs) <- 1 - transProb * (nHiddenStates - 1)
    emissionProbs <- sapply(mu, nb$density, allCounts)
    emissionProbs <- apply(emissionProbs, 2, function(v) v / sum(v)) # HMM::viterbi expects state probabilities to sum to 1 # nolint
    list(
        data = as.character(data), # HMM::viterbi expects everything codifed as strings
        HMM = HMM::initHMM(
            States = as.character(1:nHiddenStates),
            Symbols = as.character(allCounts),
            startProbs = startProbs,
            transProbs = transProbs,
            emissionProbs = t(emissionProbs)
        )
    )
}

# single-step call to initialize and run HMM::viterbi
# see initHMM for argument descriptions
viterbi.nbinomCounts <- function(nb, data, startProbs = NULL, transProb = 1e-6,
                                 mu = NULL, modeledCN = NULL, minCN = 0, maxCN = 6,
                                 asRle = TRUE # report results as runs of contiguous values using rle
                                 
){
    hmm <- initHMM(nb, data, startProbs, transProb, mu, modeledCN, minCN, maxCN)
    stateIs <- as.integer( HMM::viterbi(hmm$HMM, hmm$data) )
    results <- if(is.null(mu)) (minCN:maxCN)[stateIs] else stateIs
    if(asRle) rle(results) else results
}
