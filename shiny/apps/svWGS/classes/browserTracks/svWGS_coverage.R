#----------------------------------------------------------------------
# handle svWGS window coverage loading
#----------------------------------------------------------------------
svWGS_loadSourceCoverage <- function(sourceId){
    sessionCache$get(
        'svWGS_Coverage', 
        key = sourceId, 
        permanent = FALSE, # already an RDS file 
        from = "ram",
        create = "asNeeded",
        createFn = function(...) {
            readRDS(getSourceFilePath(sourceId, "coverageFile"))
        }  
    )$value
}
svWGS_loadSampleCoverage <- function(sourceId, sample){
    req(sourceId)
    startSpinner(session, message = "loading coverage")
    sessionCache$get(
        'svWGS_Coverage', 
        keyObject = list(sourceId, sample), 
        permanent = TRUE, 
        from = "ram",
        create = "asNeeded",
        createFn = function(...) {
            coverage <- svWGS_loadSourceCoverage(sourceId)
            startSpinner(session, message = "rebinning coverage")
            sampleCoverage <- coverage[[sample]]
            chroms <- readRDS(getSourceFilePath(sourceId, "chromosomesFile"))
            x0 <- coverage[, .(
                chrom = chrom,
                start = start + 1,
                genomeStart = getSignedNode(chroms$chromSizes, unlist(chroms$chromIndex[chrom]), start, "+", 1),
                coverage = ifelse(excluded > 0, NA, sampleCoverage) 
            )]
            svx_rebinCoverage(x0)
        }
    )$value
}

# Classes ‘data.table’ and 'data.frame':  44983 obs. of  6 variables:
#  $ chrom      : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ start      : int  0 65536 131072 196608 262144 327680 393216 458752 524288 58
# 9824 ...
#  $ gc         : num  0.449 0.378 0.499 0.451 0.393 0.477 0.392 0.429 0.468 0.462
#  ...
#  $ excluded   : int  65536 65536 65536 65536 65536 65536 65536 65536 65536 65536
#  ...
#  $ genmap     : num  0.306 0.248 0.171 0.21 0.19 ...
#  $ NA12878_cpu: num  0 0 0 0 0 0 0 0 0 0 ...
