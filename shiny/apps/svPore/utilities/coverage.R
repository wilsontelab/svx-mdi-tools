#----------------------------------------------------------------------
# handle window coverage loading
#----------------------------------------------------------------------
loadSampleCoverage <- function(sourceId, sample){
    req(sourceId)
    startSpinner(session, message = "loading coverage")
    svPoreCache$get(
        'svCoverage', 
        keyObject = list(sourceId, sample), 
        permanent = TRUE, 
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "loading coverage from disk")
            sampleCoverageZipFile <- expandSourceFilePath(sourceId, paste0(sample, ".coverage.zip")) # a direct copy of the input bgz file
            if(!file.exists(sampleCoverageZipFile)){
                sourceCoverageFile <- getSourceFilePath(sourceId, "coverageFile")
                dir <- basename(sourceCoverageFile)
                suppressWarnings(unzip(
                    sourceCoverageFile, 
                    overwrite = FALSE,
                    junkpaths = TRUE,
                    exdir = dir
                ))
                bgzFile <- list.files(
                    path = dir, 
                    pattern = paste0("^", sample, "\\..+\\.extract.windowCoverage.txt.bgz$"), 
                    full.names = TRUE
                )
                file.copy(
                    bgzFile, 
                    sampleCoverageZipFile, 
                    overwrite = FALSE
                )                
            }
            req(file.exists(sampleCoverageZipFile))
            x <- as.data.table(read.table( # fread balks at reading bgz directly, even when renamed
                sampleCoverageZipFile,
                col.names = c("chrom","start","index","coverage"),
                colClasses = c("character","integer","integer64","integer")
            ))
            binSize <- x[1:2, diff(start)]
            x[, genomeStart := index * binSize]
            x
        }
    )$value
}

# filter coverage profile by browser coordinates
filterCoverageByRange <- function(sourceId, sample, coord){
    coverage <- loadSampleCoverage(sourceId, sample)   
    startSpinner(session, message = "filtering coverage to window")
    isWholeGenome <- coord$chromosome == "all" 
    x <- if(isWholeGenome) coverage[between(as.numeric(genomeStart), coord$start, coord$end)]
                      else coverage[between(start,                   coord$start, coord$end) & chrom == coord$chromosome]
    x[, i := 1:.N]
    x
}

# aggregate coverage bins to a maximum number of points
aggregateCoverageBins <- function(coverage, maxPoints){
    coverage[, 
        point := floor(i / (ceiling(.N / maxPoints)))
    ][, 
        .(
            coverage = mean(coverage),
            index = min(index),
            genomeStart = min(genomeStart)
        ), 
        by = .(point)
    ]
}
