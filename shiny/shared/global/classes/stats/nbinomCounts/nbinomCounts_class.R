#----------------------------------------------------------------------
# object class for fitting count data to a negative binomial distribution
# an excellent model for an overdispersed Poisson distribution
#----------------------------------------------------------------------
# as used in R, the variance of the NBD is: mu + mu ** 2 / theta
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN CLASS
#----------------------------------------------------------------------
new_nbinomCounts <- function(counts) { 
    class <- 'nbinomCounts' # for reportProgress tracing

#----------------------------------------------------------------------
# execute the negative binomial fit
#----------------------------------------------------------------------
nonNegInt  <- function(counts) { # coerce to non-negative integers
    counts[is.na(counts)] <- 0
    counts[counts < 0] <- 0 
    as.integer( round(counts, 0) )
}
nonNegInt_ <- nonNegInt(counts)
allCounts <- 0:max(nonNegInt_)
nb <- MASS::glm.nb(nonNegInt_ ~ 1) # fit a negative binomial to the count data, using a single predictor class
theta <- nb$theta # shape/size/overdisperion factor, reciprocal of alpha
mu_ <- exp(coef(nb))
density <- function(mu=NULL, counts=NULL) { # calculate the probability density based on varying values of mu
    if(is.null(mu)) mu <- mu_
    if(is.null(counts)) counts <- allCounts
    dnbinom(counts, size = theta, mu = mu)
}

#----------------------------------------------------------------------
# set the return value
#----------------------------------------------------------------------
rm(nonNegInt_, nb)
structure(
    list(
        mu = mu_, # similar to Poisson lambda, suitable for dnbinom  
        theta = theta, # the shape parameter of the negative binomial, called size in dnbinom
        peak = allCounts[ which.max( density() ) ], # count at the peak of the probability density distribution
        nonNegInt = nonNegInt, # see functions above
        density = density
    ),
    class = class
)

#----------------------------------------------------------------------
# END CLASS
#----------------------------------------------------------------------
}
#----------------------------------------------------------------------
