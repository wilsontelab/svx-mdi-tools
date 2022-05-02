#----------------------------------------------------------------------
# object class for fitting count data to a negative binomial distribution
# with GC bias correction via glm.nb fit of bin read count to fraction GC
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
new_nbinomCountsGC <- function(binCounts,  # y = vector of bin read counts (class will force to non-negative integers)
                               fractionGC, # x = the fraction GC (not percent) of each bin with a count
                               binCN = 2,  # a vector of known copy numbers to stratify fit; bins with CN==NA or CN<=0 are not used during fitting # nolint
                               nGcSteps = 200, # i.e., default to 1/2 GC% steps
                               method = c('cubic', 'quadratic'), # the type of line we fit to create the model
                               maxit = 100 # value passed to glm.nb:glm.control
){ 
    class <- 'nbinomCountsGC' # for reportProgress tracing
    
#----------------------------------------------------------------------
# parse input values
#----------------------------------------------------------------------
binCounts[binCounts < 0] <- NA
binCN[binCN <= 0] <- NA # NA masks bins with missing values
readsPerAllele <- as.integer( round(binCounts / binCN, 0) )
allValues <- 0:max(readsPerAllele, na.rm = TRUE) # all allelic read counts within the range of the input data
indexGC <- round(fractionGC * nGcSteps, 0) # integer GC values of each input bin
minIndexGC <- min(indexGC, na.rm = TRUE)
maxIndexGC <- max(indexGC, na.rm = TRUE)
allIndexGC <- minIndexGC:maxIndexGC # all GC steps from min to max value
fractionGC2 <- fractionGC ** 2
fractionGC3 <- fractionGC ** 3
fractionGC_ <- allIndexGC / nGcSteps

#----------------------------------------------------------------------
# execute the negative binomial fit
#----------------------------------------------------------------------
formula <- if(method[1] == 'cubic') readsPerAllele ~ fractionGC + fractionGC2 + fractionGC3
                               else readsPerAllele ~ fractionGC + fractionGC2     
nb <- tryCatch(
    { suppressWarnings(MASS::glm.nb(formula, control = glm.control(maxit = maxit))) },
    error = function(e) {
        x <- glm(formula, family = poisson, control = glm.control(maxit = maxit))
        x$converged <- FALSE
        x$theta <- 1000
        x
    }
)
converged <- nb$converged
theta <- nb$theta # shape/size/overdisperion factor, reciprocal of alpha

#----------------------------------------------------------------------
# assemble a model summary table
#----------------------------------------------------------------------
model <- data.table(
    fractionGC = fractionGC_,    
    readsPerAllele_mu = sapply(fractionGC_, function(fgc){ # the fit value of mu as a function of binned fraction GC
        exp(predict(nb, data.frame( # exp accounts for log link in glm.nb
            fractionGC  = fgc,
            fractionGC2 = fgc ** 2,
            fractionGC3 = fgc ** 3
        )))
    })
)
model[, readsPerAllele_peak := sapply(readsPerAllele_mu, function(mu){ # the peak of the distribution for nicer visual presentation # nolint
    density <- dnbinom(allValues, size = theta, mu = mu) 
    allValues[ which.max( density ) ]
})]

#----------------------------------------------------------------------
# set the return value
#----------------------------------------------------------------------
structure(
    list(
        nGcSteps = nGcSteps,
        converged = converged,
        theta = theta, # the shape parameter of the negative binomial, called size in dnbinom        
        model = model, # model data, one row per GC bin with mu and peak readsPerAllele values
        indexGC_offset = minIndexGC - 1
    ),
    class = class
)

#----------------------------------------------------------------------
# END CLASS
#----------------------------------------------------------------------
}
#----------------------------------------------------------------------
