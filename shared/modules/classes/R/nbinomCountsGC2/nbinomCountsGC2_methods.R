#----------------------------------------------------------------------
# nbinomCountsGC2 class generic methods, called as method(obj)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# returns readsPerAllele from a set of input GC values (or the model itself)
#----------------------------------------------------------------------
getRows.nbinomCountsGC2 <- function(nb, fractionGC){
    if(is.null(fractionGC)) TRUE else {
        pmax(nb$minGcIndex, pmin(nb$maxGcIndex, round(fractionGC * nb$nGcSteps, 0))) - nb$gcIndexOffset
    }
}
predict.nbinomCountsGC2 <- function(
    nb, 
    fractionGC = NULL, 
    type = c('mu', 'peak', 'adjustedPeak', 'theta'), 
    peakThreshold = 10
){
    rows <- getRows.nbinomCountsGC2(nb, fractionGC)
    if(type[1] == 'adjustedPeak'){ # peak (mode) is best for visualization
        peak <- nb$peak[rows]  
        mu   <- nb$mu[rows]
        ifelse(peak < peakThreshold, mu, peak) # revert to mu when peak is unreliable due to collision with zero counts
    } else {  
        nb[[ type[1] ]][rows] # default is mu, best for HMM and similar
    }
}

#----------------------------------------------------------------------
# returns a vector of cumulative count probabilities for a set of bins relative to a computed model
# each bin's value is its lower-tail probability, i.e., percentile of its actual count based on the nb model
#----------------------------------------------------------------------
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

#----------------------------------------------------------------------
# solve for the most likely CN path of a set of bin counts given a nbinomCountsGC2 model
# depends on class hmmEPTable
#----------------------------------------------------------------------
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
    rows <- getRows.nbinomCountsGC2(nb, fractionGC)
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

#----------------------------------------------------------------------
# make a composite plot of the GC bias model for cell quality monitoring
#----------------------------------------------------------------------
plot.nbinomCountsGC2 <- function(nb, gc_w, NR_map_w, modal_CN = 2, rejected = FALSE, 
                                 col = NULL, maxPoints = 5000, binCN = NULL){
    peak <- predict(nb, type = 'adjustedPeak') * modal_CN
    maxPeak  <- max(peak, na.rm = TRUE)
    if(!is.null(binCN)) NR_map_w <- NR_map_w / binCN * modal_CN
    maxCount <- max(NR_map_w, na.rm = TRUE)
    ymax <- min(maxPeak * 2, maxCount)
    n <- length(gc_w)
    i <- sample(n, min(n, maxPoints), replace = FALSE)
    col <- if(is.null(col) || is.na(col)) rgb(0, 0, 0, 0.1) else col[i]
    plot(gc_w[i], NR_map_w[i], 
         pch = 19, cex = 0.4, col = col,
         xlim = c(0.3, 0.6), ylim = c(0, ymax),
         xlab = "Fraction GC", ylab = "# of Reads")
    col <- if(rejected) "red3" else "blue"
    lines(
        nb$gcFractions,
        peak,
        lty = 1, lwd = 1.5, col = col
    )
    for(percentile in c(0.05, 0.95)) lines(
        nb$gcFractions, 
        qnbinom(percentile, size = nb$theta, mu = nb$mu) * modal_CN, 
        lty = 3, lwd = 1.5, col = col
    )
}
