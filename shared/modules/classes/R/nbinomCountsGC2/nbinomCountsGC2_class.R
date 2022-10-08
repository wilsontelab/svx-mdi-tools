#----------------------------------------------------------------------
# object class for fitting count data to a negative binomial distribution
# with GC bias correction via glm.nb fit of bin read count to fraction GC
#----------------------------------------------------------------------
# in version 2 of the model, i.e., nbinomCountsGC2, both mu and theta
# are allowed to vary as a function of GC + GC**2
#----------------------------------------------------------------------
# as used in R, the variance of the NBD is: mu + mu**2 / theta
#----------------------------------------------------------------------
# read counts are rounded to integers as needed
# negative read counts are set to NA
# NA read counts are ignored in model building by glm.nb
# if requested, only some bins with known copy number can be used to fit the model
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN CLASS
#----------------------------------------------------------------------
new_nbinomCountsGC2 <- function(
    binCounts,  # y = vector of bin read counts (class will force to non-negative integers)
    fractionGC, # x = the fraction GC (not percent) of each bin with a count
    binCN = 2,  # a vector of known copy numbers to stratify fit; bins with CN==NA or CN<=0 are not used during fitting # nolint
    nGcSteps = 100, # i.e., default to 1 GC% steps
    maxit = 100 # value passed to glm.nb:glm.control
){ 
    class <- 'nbinomCountsGC2' # for reportProgress tracing
    
#----------------------------------------------------------------------
# parse input values
#----------------------------------------------------------------------
binCounts[binCounts < 0] <- NA
binCN[binCN <= 0] <- NA # NA masks bins with missing values
readsPerAllele <- as.integer( round(binCounts / binCN, 0) )
allValues <- 0:max(readsPerAllele, na.rm = TRUE) # all allelic read counts within the range of the input data
gcIndex <- round(fractionGC * nGcSteps, 0) # integer GC values of each input bin
minGcIndex <- min(gcIndex, na.rm = TRUE)
maxGcIndex <- max(gcIndex, na.rm = TRUE)
gcIndices  <- minGcIndex:maxGcIndex # all GC steps in use from min to max value
gcIndices2 <- gcIndices ** 2

#----------------------------------------------------------------------
# execute the negative binomial fits across all GC bins
#----------------------------------------------------------------------
fits <- sapply(gcIndices, function(gci){
    rpa <- readsPerAllele[gcIndex == gci]
    rpa <- rpa[!is.na(rpa)]
    n <- length(rpa)
    if(n < 5) return(c(1, 10, 0)) # unweight rare bins on edges of gc space
    tryCatch(
        { 
            x <- suppressWarnings(MASS::glm.nb(rpa ~ 1, control = glm.control(maxit = maxit))) 
            c(exp(coef(x)), x$theta, n)
        },
        error = function(e) {
            x <- glm(rpa ~ 1, family = "poisson", control = glm.control(maxit = maxit))
            c(exp(coef(x)), 10000, n) # assume high theta, i.e, not over-dispersed
        }
    )
})

#----------------------------------------------------------------------
# calculate the fitted mu and theta values for all GC index values
#----------------------------------------------------------------------
mu <- fits[1, ]
mu <- lm(mu ~ gcIndices + gcIndices2, weights = fits[3, ])
mu <- predict(mu)
mu <- ifelse(mu < 1, median(mu, na.rm =TRUE), mu)
mu <- pmax(0.1, mu)
theta <- pmax(0.1, fits[2, ])
peak <- mapply(
    function(size, mu) tryCatch({
        density <- dnbinom(allValues, size = size, mu = mu)
        allValues[ which.max( density ) ]
    }, error = function(e){
        message(size)
        message(mu)
        stop("bad value for mu or theta: nbinomCountsGC2 peak")
    }), 
    size = theta, 
    mu = mu
)

#----------------------------------------------------------------------
# set the return value
#----------------------------------------------------------------------
structure(
    list(
        nGcSteps = nGcSteps, # value for parsing new GC values to mu and theta values
        gcIndexOffset = minGcIndex - 1,
        gcFractions = gcIndices / nGcSteps,
        theta = theta, # the shape parameter of the negative binomial, called size in dnbinom             
        mu = mu,       # poisson mean, as reads per allele, used by HMM
        peak = peak    # the most frequent reads per allele value in each GC bin, best for visualization
    ),
    class = class
)

#----------------------------------------------------------------------
# END CLASS
#----------------------------------------------------------------------
}
#----------------------------------------------------------------------
