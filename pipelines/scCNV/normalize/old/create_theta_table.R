#=====================================================================================
# this script is not called by the pipeline
# it was run once to generate the lookup table for estimating theta ~ ER_modal_CN + sdLagDiff
# in other words, to get an over-dispersion estimate based on the inter-window lagging difference
# scaled to the peak window value as an estimate of mu
#=====================================================================================

# set a range of working mu values, ER_modal_CN is matched to one of them
mus <- data.frame(
    mu = as.vector(sapply(1:5, function(exp) 10**exp * seq(1, 9.5, 0.5))),
    start = NA,
    end = NA
)
mus$end   <- mus$mu + c(diff(mus$mu) * 0.5, Inf)
mus$start <- c(0, mus$end[1:(nrow(mus) - 1)])

# set a range of working theta values, ultimately used for generating count likelihoods
thetas <- 10**seq(-3, 5, 0.5)

# use random simulations to estimate sdLagDiff for every combination of mu and theta
nIter <- 100000
sdLagDiffs <- sapply(thetas, function(theta){ # yields a matrix[mu, theta]
    sapply(df$mu, function(mu){
        counts <- rnbinom(nIter, size = theta, mu = mu)
        sd(diff(counts))
    })
})

# create a function for doing the lookup
getTheta <- function(x, ER_modal_CN, sdLagDiff) with(x, {
    i <- which(ER_modal_CN >= mus$start & ER_modal_CN < mus$end)
    j <- which.min(abs(sdLagDiff - sdLagDiffs[i, ]))
    thetas[j]
})

# save all required objects for loading and theta estimation
# usage:
#    thetaLookup <- readRDS(path)
#    theta <- thetaLookup$getTheta(thetaLookup, ER_modal_CN, sdLagDiff)
thetaLookup <- list(
    mus = mus,
    thetas = thetas,
    sdLagDiffs = sdLagDiffs,
    getTheta = getTheta
)

dir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/pipelines/scCNV/normalize"
file <- file.path(dir, "theta_lookup_table.rds")
saveRDS(thetaLookup, file = file)
