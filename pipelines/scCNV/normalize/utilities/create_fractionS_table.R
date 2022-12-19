#=====================================================================================
# this script is not called by the pipeline
# it was run once to generate a lookup table for establishing fraction of mono-allelic
# replication states (i.e., N = 3) as a function of fractionS
#=====================================================================================

# model parameters
nLoci <- 10000
weigthPower <- 5 # how much to prefer re-replication, rather than replicating a new window
                 # larger numbers assume increasingly synchronous execution of two alleles when ploidy==2
                 # if anything, could still move from 5 down to 4 for less synchrony

# randomly select allelic replication order
weights <- (nLoci:1) ** weigthPower
x <- matrix(NA, nrow = 2, ncol = nLoci)
x[1, ] <- sample(nLoci, prob = weights) # interleave weighted window choices
x[2, ] <- sample(nLoci, prob = weights)
x <- as.vector(x)

# assemble the table correlating fractionS to window replication states
fractionS <- seq(0.05, 0.95, 0.025)
table <- data.frame(
    fractionS = fractionS,
    fractionS_int = as.integer(fractionS * 1000),
    P1_R0 = 1 - fractionS,
    P1_R1 = fractionS
)
table <- cbind(table, t(sapply(table$fractionS, function(fractionS){
    xx <- x[1:(nLoci * fractionS * 2)]
    agg <- aggregate(xx, list(xx), length)
    agg <- aggregate(agg[[2]], list(agg[[2]]), length)
    n1 <- agg[[1]] == 1
    n1 <- if(any(n1)) agg[n1, 2] else 0
    n2 <- agg[[1]] == 2
    n2 <- if(any(n2)) agg[n2, 2] else 0
    c(nLoci - n1 - n2, n1, n2) / nLoci
})))
names(table) <- c("fractionS", "fractionS_int", "P1_R0", "P1_R1", "P2_R0", "P2_R1", "P2_R2")
table <- rbind(
    data.frame(
        fractionS = 0,
        fractionS_int = 0L,
        P1_R0 = 1,
        P1_R1 = 0,
        P2_R0 = 1,
        P2_R1 = 0,
        P2_R2 = 0
    ),
    table
)

# create a function for doing the lookup
getRepStateProbs <- function(x, ploidy, fractionS) with(x, {
    fractionS_int <- as.integer(fractionS * 1000)
    row <- table[table$fractionS_int == fractionS_int, , drop = FALSE]    
    # if(nrow(row) > 0){ # a perfect match
        if(ploidy == 1) c(row$P1_R0,         0, row$P1_R1) 
                   else c(row$P2_R0, row$P2_R1, row$P2_R2)        
    # } else { # in between value, interpolate
    #     below <- table[max(which(table$fractionS_int < fractionS_int)), , drop = FALSE]    
    #     above <- table[min(which(table$fractionS_int > fractionS_int)), , drop = FALSE]    
    #     if(ploidy == 1) c(mean(below$P1_R0, above$P1_R0),                              0, mean(below$P1_R1, above$P1_R1)) 
    #                else c(mean(below$P2_R0, above$P2_R0), mean(below$P2_R1, above$P2_R1), mean(below$P2_R2, above$P2_R2))
    # }
})

# save the table for use during fitting
# usage:
#    fractionSLookup <- readRDS(path)
#    probs <- fractionSLookup$getRepStateProbs(fractionSLookup, ploidy, fracS) # for 0, 1, and 2 alleles replicated
fractionSLookup <- list(
    table = table,
    getRepStateProbs = getRepStateProbs
)
dir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/mdi/suites/definitive/svx-mdi-tools/pipelines/scCNV/normalize"
file <- file.path(dir, "fractionS_table.rds")
saveRDS(fractionSLookup, file = file)
