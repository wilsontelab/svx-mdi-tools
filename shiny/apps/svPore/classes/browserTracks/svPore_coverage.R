#----------------------------------------------------------------------
# handle svPore window coverage loading
#----------------------------------------------------------------------
svPore_loadSampleCoverage <- function(sourceId, sample){
    req(sourceId)
    startSpinner(session, message = "loading coverage")
    sessionCache$get(
        'svPore_Coverage', 
        keyObject = list(sourceId, sample), 
        permanent = TRUE, 
        from = "ram",
        create = "asNeeded",
        createFn = function(...) {
            startSpinner(session, message = "unpacking coverage .")
            sampleCoverageZipFile <- expandSourceFilePath(sourceId, paste0(sample, ".coverage.zip")) # a direct copy of the input bgz file
            if(!file.exists(sampleCoverageZipFile)){
                sourceCoverageFile <- getSourceFilePath(sourceId, "coverageFile")
                dir <- dirname(sourceCoverageFile)
                startSpinner(session, message = "unpacking coverage ..")
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
                startSpinner(session, message = "unpacking coverage ...")
                file.copy(
                    bgzFile, 
                    sampleCoverageZipFile, 
                    overwrite = FALSE
                )                
            }

            startSpinner(session, message = "rebinning coverage")
            req(file.exists(sampleCoverageZipFile))
            # TODO: convert this to bin coverage, like svWGS?
            x0 <- as.data.table(read.table( # fread balks at reading bgz directly, even when renamed
                sampleCoverageZipFile,
                col.names = c("chrom","start","index","coverage"),
                colClasses = c("character","integer","integer64","integer")
            ))
            binSize <- x0[1:2, diff(start)]
            x0[, ":="(
                genomeStart = index * binSize,
                gc = NA_real_,
                excluded = 0
            )]
            svx_rebinCoverage(sourceId, sample, x0)
        }
    )$value
}
