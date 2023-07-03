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
                dir <- dirname(sourceCoverageFile)
                startSpinner(session, message = "unzipping from package")
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
            startSpinner(session, message = "reading from disk")
            req(file.exists(sampleCoverageZipFile))
            x0 <- as.data.table(read.table( # fread balks at reading bgz directly, even when renamed
                sampleCoverageZipFile,
                col.names = c("chrom","start","index","coverage"),
                colClasses = c("character","integer","integer64","integer")
            ))
            startSpinner(session, message = "aggregating bins")
            binSize <- x0[1:2, diff(start)]
            binDensity <- 1 / binSize
            medianCoverage <- x0[, median(coverage)]
            aggFactor <- 11
            x0[, ":="(
                genomeStart = index * binSize,
                index1 = floor(index / aggFactor) # for the next round aggregation
            )]
            x1 <- x0[, .(
                start = start[5], # start of the middle bin of an 11-bin span
                genomeStart = genomeStart[5],            
                coverage = mean(coverage),
                index2 = floor(index1 / aggFactor)
            ), by = .(chrom, index1)]
            x2 <- x1[, .(
                start = start[5],
                genomeStart = genomeStart[5],            
                coverage = mean(coverage)
            ), by = .(chrom, index2)]
            outCols <- c("chrom","start","genomeStart","coverage")
            list(
                binSize = binSize,
                binDensities = sapply(0:2, function(power) binDensity / aggFactor**power),
                medianCoverage = medianCoverage,
                x0 = x0[, .SD, .SDcols = outCols],
                x1 = x1[, .SD, .SDcols = outCols],
                x2 = x2[, .SD, .SDcols = outCols]
            )
        }
    )$value
}

# filter coverage profile by browser coordinates
filterCoverageByRange <- function(sourceId, sample, coord, maxPoints){
    coverage <- loadSampleCoverage(sourceId, sample)   
    startSpinner(session, message = "filtering coverage to window")
    maxBinDensity <- maxPoints / coord$width
    passingPowers <- which(coverage$binDensities <= maxBinDensity) - 1L
    power <- if(length(passingPowers) == 0) 2
        else if(passingPowers[1] == 0) 0 
        else passingPowers[1] - 1 # thus, if possible, return the first density that did NOT pass, to allow partial aggregation below
    x <- coverage[[paste0("x", power)]]
    isWholeGenome <- coord$chromosome == "all" 
    x <- if(isWholeGenome) x[between(as.numeric(genomeStart), coord$start, coord$end)]
                      else x[between(start,                   coord$start, coord$end) & chrom == coord$chromosome]
    x[, i := 1:.N]
    list(
        binSize = coverage$binSize,
        power = power,
        medianCoverage = coverage$medianCoverage,
        x = x
    )
}

# aggregate coverage bins to a maximum number of points
aggregateCoverageBins <- function(coverage, maxPoints){
    coverage
    startSpinner(session, message = "aggregating bins")
    coverage$x[, point := floor(i / (ceiling(.N / maxPoints)))]
    coverage$x <- coverage$x[, 
        .(
            start = median(start),
            genomeStart = median(genomeStart),           
            coverage = mean(coverage)
        ), 
        by = .(chrom, point)
    ]
    coverage
}

# support read depth and copy number plots
setCoverageXY <- function(coverage_, plotAs, medianPloidy, coord, sample){
    switch(
        plotAs,
        "Read Depth" = coverage_$x[, .(
            x = if(coord$chromosome == "all") as.numeric(genomeStart) else start,
            y = coverage,
            sample = sample
        )],
        coverage_$x[, .( # Copy Number (possibly converted later to Copy Number Change)
            x = if(coord$chromosome == "all") as.numeric(genomeStart) else start,
            y = coverage / coverage_$medianCoverage * medianPloidy,
            sample = sample
        )]
    )
}
