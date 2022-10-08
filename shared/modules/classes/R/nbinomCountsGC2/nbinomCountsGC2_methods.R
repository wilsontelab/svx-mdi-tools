#----------------------------------------------------------------------
# nbinomCountsGC2 class generic methods, called as method(obj)
#----------------------------------------------------------------------

# returns readsPerAllele from a set of input GC values (or the model itself)
predict.nbinomCountsGC2 <- function(
    nb, 
    fractionGC = NULL, 
    type = c('mu', 'peak', 'adjustedPeak', 'theta'), 
    peakThreshold = 10
){
    rows <- if(is.null(fractionGC)) TRUE else round(fractionGC * nb$nGcSteps, 0) - nb$gcIndexOffset    
    if(type[1] == 'adjustedPeak'){
        peak <- nb$peak[rows]  
        mu   <- nb$mu[rows]
        ifelse(peak < peakThreshold, mu, peak) # revert to mu when peak is unreliable due to collision with zero counts
    } else {  
        nb[[ type[1] ]][rows]        
    }
}

# returns a vector of cumulative count probabilities for a set of bins relative to a computed model
# each bin's value is its lower-tail probability, i.e., percentile for its actual count based on the nb model
cumprob <- function(x, ...) {
    UseMethod("cumprob", x)
}
cumprob.nbinomCountsGC2 <- function(
    nb,
    binCounts,  # vector of bin read counts (will force to non-negative integers)
    fractionGC, # the fraction GC (not percent) of each bin with a count
    binCN = 2   # a vector of known copy numbers to stratify fit; bins with CN==NA or CN<=0 are not used during fitting # nolint
){
    binCounts[binCounts < 0] <- NA
    binCN[binCN <= 0] <- NA # NA masks bins with missing values
    readsPerAllele <- as.integer( round(binCounts / binCN, 0) )
    allValues <- 0:max(readsPerAllele, na.rm = TRUE)    
    gcIndex <- round(fractionGC * nb$nGcSteps, 0) - nb$gcIndexOffset
    lowertail <- list() # buffer to avoid re-calling dnbinom
    sapply(seq_along(readsPerAllele), function(i){
        if(is.na(readsPerAllele[i])) return(NA)
        j <- gcIndex[i]
        if(is.null(unlist(lowertail[j]))) lowertail[[j]] <<- {
            pr <- dnbinom(allValues, size = nb$theta[j], mu = nb$mu[j])
            cumsum(pr)            
        }
        lowertail[[j]][readsPerAllele[i] + 1]
    })
}

# solve for the most likely CN path of a set of bin counts given a nbinomCountsGC2 model
viterbi <- function(x, ...) {
    UseMethod("viterbi", x)
}
viterbi.nbinomCountsGC2 <- function(
    nb, # a nbinomCountsGC2 model
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
    rows <- round(fractionGC * nb$nGcSteps, 0) - nb$gcIndexOffset
    size <- nb$theta[rows]
    rpa  <- nb$mu[rows]
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

# save a plot of the GC bias correction for quality monitoring
plot.nbinomCountsGC2 <- function(nb, cell, ploidy, filename){
    png(
        filename = filename,
        width = 1.5, 
        height = 1.5, 
        units = "in", 
        pointsize = 6,
        bg = "white",  
        res = 300,
        type = "cairo"
    )
    par(mar= c(4.1, 4.1, 0.1, 0.1))
    plot(cell$gc_wr, cell$NR_map_wr, xlim = c(0.35, 0.55), 
         pch = 16, cex = 0.25, col = rgb(0, 0, 0, 0.2))
    col <- if(cell$rejected) "red" else "blue"
    lines(
        nb$gcFractions,
        predict(nb, type = 'adjustedPeak') * ploidy,
        lty = 1, lwd = 1, col = col
    )
    for(percentile in c(0.05, 0.95)) lines(
        nb$gcFractions, 
        qnbinom(percentile, size = nb$theta, mu = nb$mu) * ploidy, 
        lty = 3, lwd = 0.75, col = col
    )
    dev.off()
}
