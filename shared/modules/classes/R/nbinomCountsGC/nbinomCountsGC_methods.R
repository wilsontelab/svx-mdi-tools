#----------------------------------------------------------------------
# nbinomCountsGC class generic methods, called as method(obj)
#----------------------------------------------------------------------

# returns an object with both x and y data (since x may have been our default) binned fraction GC
predict.nbinomCountsGC <- function(
    nb, 
    fractionGC=NULL, 
    type=c('mu', 'peak', 'adjustedPeak'), 
    peakThreshold = 10
){
    rows <- if(is.null(fractionGC)) TRUE else round(fractionGC * nb$nGcSteps, 0) - nb$indexGC_offset    
    if(type[1] == 'adjustedPeak'){
        peak <- nb$model$readsPerAllele_peak[rows]  
        mu   <- nb$model$readsPerAllele_mu[rows]
        ifelse(peak < peakThreshold, mu, peak) # revert to mu when peak is unreliable due to collision with zero
    } else {
        col <- paste('readsPerAllele', type[1], sep = '_')    
        nb$model[[col]][rows]        
    }
}

# returns a vector of cumulative count probabilities for a set of bins relative to a computed model
# each bin's value is its lower-tail probability for its actual count based on the nb model for its GC fraction
cumprob <- function(x, ...) {
    UseMethod("cumprob", x)
}
cumprob.nbinomCountsGC <- function(
    nb,
    binCounts,  # vector of bin read counts (will force to non-negative integers)
    fractionGC, # the fraction GC (not percent) of each bin with a count
    binCN = 2   # a vector of known copy numbers to stratify fit; bins with CN==NA or CN<=0 are not used during fitting # nolint
){
    binCounts[binCounts < 0] <- NA
    binCN[binCN <= 0] <- NA # NA masks bins with missing values
    readsPerAllele <- as.integer( round(binCounts / binCN, 0) )
    allValues <- 0:max(readsPerAllele, na.rm = TRUE)    
    indexGC <- round(fractionGC * nb$nGcSteps, 0) - nb$indexGC_offset
    lowertail <- list() # buffer to avoid re-calling dnbinom
    sapply(seq_along(readsPerAllele), function(i){
        if(is.na(readsPerAllele[i])) return(NA)
        j <- indexGC[i]
        if(is.null(unlist(lowertail[j]))) lowertail[[j]] <<- {
            pr <- dnbinom(allValues, size = nb$theta, mu = nb$model$readsPerAllele_mu[j])
            cumsum(pr)            
        }
        lowertail[[j]][readsPerAllele[i] + 1]
    })
}

# solve for the most likely CN path of a set of bin counts given a nbinomCountsGC model
viterbi <- function(x, ...) {
    UseMethod("viterbi", x)
}
viterbi.nbinomCountsGC <- function(
    nb, # a nbinomCountsGC model
    binCounts, 
    fractionGC, # similar to original arguments of constructor
    percentile = NULL, # bin medians, optional result of a batch effect normalization, to adjust mu=readsPerAllele # nolint
    maxCN = 6, 
    transProb = 1e-6, # options for the HMM
    chroms = NULL, # if a vector, use keyedViterbi by chromosome
    forceCNs = NULL, # if a vector, force to HMM output for bins that are not NA
    asRle = TRUE
){ # report results as rle object to minimize object size

    # calculate emissProbs
    size <- nb$theta 
    rpa <- predict(nb, fractionGC, type = 'mu')
    if(!is.null(percentile)) rpa <- qnbinom(percentile, size = size, mu = rpa)     
    CNs <- 0:maxCN    
    binCounts[binCounts < 0] <- NA
    binCounts <- as.integer( round(binCounts, 0) )  
    emissProbs <- sapply(CNs, function(CN){
        dnbinom(binCounts, size = size, mu = if(CN == 0) 0.1 else rpa * CN, log = TRUE)
    })

    # apply forced states
    forceStateIs <- if(!is.null(forceCNs)) {
        forceCNs <- pmin(maxCN, forceCNs)
        forceCNs[forceCNs < 0] <- NA
        forceCNs + 1
    } else NULL

    # construct and solve the hmm
    hmm <- new_hmmEPTable(
        emissProbs,
        transProb,
        keys = chroms,
        forceStateIs = forceStateIs
    )
    stateIs <- keyedViterbi(hmm) # reverts to viterbi over all observations if chroms==NULL

    # parse and return the results as CN
    results <- CNs[stateIs]
    list(
        hmm = hmm,
        cn = if(asRle) rle(results) else results,
        maxCN = maxCN
    )
}
