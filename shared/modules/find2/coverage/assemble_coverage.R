# merge samples and add normalization data

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(yaml)
    library(bit64)
    library(parallel)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOMES_DIR',
        'GENOME',
        'FIND_PREFIX',
        'FIND_MODE',
        'TASK_DIR',
        'OUTPUT_DIR',
        'DATA_NAME'
    ), 
    integer = c(
        'N_CPU',
        'BIN_SIZE',
        'KMER_LENGTH',
        'N_ERRORS'
    )
))
env$ASSESS_CNV_BINS <- Sys.getenv("ASSESS_CNV_BINS") != ""
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(rUtilDir, c('utilities'))
sourceScripts(file.path(rUtilDir, 'genome'), c('general', 'chroms'))
#-------------------------------------------------------------------------------------
# parse the project name and directory
if(env$FIND_MODE == "find"){
    projectName <- basename(env$OUTPUT_DIR)
    projectDir  <- env$OUTPUT_DIR
    getSamplePrefix <- function(sample) file.path(env$TASK_DIR, sample)
} else {
    projectName <- env$DATA_NAME
    projectDir  <- env$TASK_DIR
    getSamplePrefix <- function(sample) file.path(env$TASK_DIR, sample, sample)
}
#=====================================================================================

#=====================================================================================
# load the metadata across all samples
#-------------------------------------------------------------------------------------
message("loading sample metadata")
inFile <- paste(env$FIND_PREFIX, "metadata", "yml", sep=".")
metadata <- read_yaml(inFile)
metadata <- lapply(metadata, function(x) strsplit(as.character(x), "\\s+")[[1]])
#=====================================================================================

#=====================================================================================
# load fixed-width bins data
#-------------------------------------------------------------------------------------
message("loading bin metadata")
binsDir <- file.path(env$GENOMES_DIR, "bins", env$GENOME, "fixed_width_bins")
chunkSize <- 65536 # fixed by extract/base_coverage.pl
binsFile <- paste0(env$GENOME, ".bins.size_", env$BIN_SIZE, ".k_", env$KMER_LENGTH, ".e_", env$N_ERRORS, ".bed.gz")
binsFile <- file.path(binsDir, binsFile)
bins <- if(file.exists(binsFile)) {
    fread(
        binsFile,
        sep = "\t",
        header = FALSE,
        colClasses = c(
            "character", 
            "integer", 
            "integer", 
            "integer", 
            "numeric", 
            "character", 
            "integer", 
            "integer", 
            "integer", 
            "numeric", 
            "numeric"
        ),
        col.names = c(
            "chrom", 
            "start", 
            "end", 
            "cumIndex", 
            "gc", 
            "strand", 
            "excluded", 
            "gap", 
            "bad", 
            "umap", 
            "genmap"
        )
    )[, .SD, .SDcols = c(
        "chrom",
        "start",
        "end",
        "gc",
        "excluded",
        "genmap"
    )]
} else {
    message("WARNING: missing bins file")
    message(paste0("    ", binsFile))
    message(paste0("    ", "proceeding with dummy values for gc, excluded, genmap"))
    setCanonicalChroms()
    chromSizes <- loadChromSizes(windowSize = env$BIN_SIZE)
    start <- (1:as.integer(nChromWindows) - 1) * env$BIN_SIZE
    chromSizes[, .(
        start = start,
        gc = 0.5,
        excluded = 0,
        genmap = 1 
    ), by = .(chrom)]
}
#=====================================================================================

#=====================================================================================
# load CNVs suitable for copy number characterization
#-------------------------------------------------------------------------------------
if(env$ASSESS_CNV_BINS){
    message("loading called CNVs")
    cnvs <- readRDS(paste(env$FIND_PREFIX, 'structural_variants', 'rds', sep = "."))
    cnvs <- cnvs[JXN_TYPE %in% c("D", "L")]
}
#=====================================================================================

#=====================================================================================
# support functions for binning coverage
#-------------------------------------------------------------------------------------
getSampleCoverageIndex <- function(samplePrefix){
    coverageIndexFile <- paste(samplePrefix, env$GENOME, "coverage.index.gz", sep = ".")
    coverageIndex <- fread(
        coverageIndexFile,
        sep = "\t",
        header = TRUE,
        colClasses = c(
            "integer", # chromIndex, 1-referenced
            "integer", # chunkIndex, 0-referenced
            "numeric", # cumNBreaks printed to breaks and counts files through (i.e., including) this chunk
            "numeric"  # coverage
        )
    )
    coverageIndex[, ":="(
        chrom = metadata$CHROMS[chromIndex],
        start = chunkIndex * chunkSize # chrom coordinate, 0-referenced (like BED)
    )]    
}
assignBins <- function(breakSpans){ # spans either on a whole chromosome, or within a specific CNV
    binAssignments <- rbind(
        if(breakSpans[, any(sameBin)]) breakSpans[sameBin == TRUE, .( # handle single-bin break spans separately for speed
            breakStart = breakStart,
            breakEnd   = breakEnd,
            coverage   = coverage,
            start      = bin1       
        )] else NULL, 
        if(breakSpans[, any(!sameBin)]) breakSpans[sameBin == FALSE, .( # thus, some breaks are assigned to two bins
            breakEnd  = breakEnd,
            coverage  = coverage,
            start     = seq(bin1, bin2, env$BIN_SIZE)
        ), by = .(breakStart)] else NULL
    )
    binAssignments[, end := start + env$BIN_SIZE]
}
#=====================================================================================

#=====================================================================================
# split CNV candidates against the bins for CN estimation after GC bias correction in app
#-------------------------------------------------------------------------------------
message("assembling coverage by bin")
cnvBins <- data.table(svId = character(), sample = character(), gc = double(), coverage = double(), overlap = integer())

# initialize each sample
for(sample in metadata$SAMPLES){
    message(paste("   ", sample))
    samplePrefix <- getSamplePrefix(sample)
    coverageIndex <- getSampleCoverageIndex(samplePrefix) 
    sampleCnvs <- if(env$ASSESS_CNV_BINS) cnvs[cnvs[[sample]] >= 2] else NULL # ignore singletons for efficiency, they aren't likely to be informativ re: CN
    offsets <- c(0, coverageIndex$cumNBreaks) 
    breaksFile <- paste(samplePrefix, env$GENOME, "coverage.breaks", sep = ".")
    countsFile <- paste(samplePrefix, env$GENOME, "coverage.counts", sep = ".")
    breaksFile <- file(breaksFile, "rb")
    countsFile <- file(countsFile, "rb")    
    binCoverage <- data.table(chrom = character(), start = integer(), end = integer(), coverage = double())
    cnvBinCoverage <- cbind(binCoverage, data.table(svId = character(), overlap = integer()))
    nullCnvBinCoverage <- cnvBinCoverage

    # work on chromosome at a time
    do.call(rbind, lapply(unique(coverageIndex$chrom), function(chrom_){ # this cannot be reliably parallelized, seemingly due to seek/read collisions?
        message(paste("     ", chrom_)) 

        # load the break spans, i.e., the contiguous segments of the genome covered by the same number of reads
        breakSpans <- coverageIndex[
            chrom == chrom_,
            {
                offset  <- offsets[.I]
                nBreaks <- cumNBreaks - offset
                seek(breaksFile, where = offset * 2, origin = "start", rw = "read")
                seek(countsFile, where = offset,     origin = "start", rw = "read")
                .(
                    breaks = readBin(breaksFile, "integer", n = nBreaks, size = 2, signed = FALSE) + start,
                    counts = readBin(countsFile, "integer", n = nBreaks, size = 1, signed = FALSE)
                )
            },
            by = .(chunkIndex)
        ]

        # pad the first and last break positions out to the ends of the chromosome
        breakSpans <- breakSpans[, .(
            breakStart = c(0, breaks),
            breakEnd   = c(breaks, coverageIndex[chrom == chrom_, max(start) + chunkSize]), # end will likely overrun the true chromosome a bit
            coverage   = c(0, counts)
        )]

        # assign the break span ends to their proper fixed width bins; most break spans are contained in a single bin, but not all
        breakSpans[, ":="(
            bin1 = floor(breakStart / env$BIN_SIZE) * env$BIN_SIZE,
            bin2 = floor(breakEnd   / env$BIN_SIZE) * env$BIN_SIZE
        )]
        breakSpans[, ":="(sameBin = bin1 == bin2)]

        # calculate genome coverage by bin for this sample+chrom
        binAssignments <- assignBins(breakSpans)   
        binCoverage <- rbind(binCoverage, binAssignments[, .(
            chrom = chrom_,
            coverage = weighted.mean(coverage, pmin(breakEnd, end) - pmax(breakStart, start))
        ), keyby = .(start, end)])

        # calculate genome coverage by bin, per CNV, for GC bias correction of CNV CN in app
        if(!env$ASSESS_CNV_BINS) next
        if(nrow(sampleCnvs) == 0) next
        sampleChromCnvs <- sampleCnvs[CHROM_1 == chrom_]
        nCnvs <- nrow(sampleChromCnvs)
        if(nCnvs == 0) next
        setkey(breakSpans, breakStart)
        cnvBinCoverage <<- rbind(cnvBinCoverage, do.call(rbind, mclapply(1:nCnvs, function(i){  
            cnv <- sampleChromCnvs[i]
            cnvStart <- cnv[, min(POS_1, POS_2) - 1]
            cnvEnd   <- cnv[, max(POS_1, POS_2)]
            bs <- breakSpans[breakStart %between% c(cnvStart, cnvEnd)]
            if(nrow(bs) == 0) return(nullCnvBinCoverage) # catch rare CNVs with no breaks
            binAssignments <- assignBins(bs) 
            binAssignments[, {
                overlaps <- pmin(breakEnd, end, cnvEnd) - pmax(breakStart, start, cnvStart)
                .(
                    chrom = chrom_,
                    coverage = weighted.mean(coverage, overlaps),
                    svId = cnv$SV_ID,
                    overlap = sum(overlaps)
                )
            }, keyby = .(start, end)]
        }, mc.cores = env$N_CPU)))
    }))
    close(breaksFile)
    close(countsFile)

    # add this sample's coverage to the bins table
    bins <- merge(
        bins,
        binCoverage[, .(
            chrom = chrom,
            start = start,
            coverage = coverage
        )],
        by = c("chrom", "start"),
        all.x = TRUE
    )
    names(bins)[ncol(bins)] <- sample

    # append this samples CNV bins
    if(env$ASSESS_CNV_BINS){
        cnvBinCoverage <- merge(
            cnvBinCoverage, 
            bins[, .(chrom, start, gc, excluded)],
            by = c("chrom", "start"),
            all.x = TRUE
        )
        cnvBins <- rbind(cnvBins, cnvBinCoverage[excluded / env$BIN_SIZE < 0.1, .(
            sample   = sample,
            gc       = paste(round(gc, 3),       collapse = ","),
            coverage = paste(round(coverage, 1), collapse = ","),
            overlap  = paste(overlap,            collapse = ",")
        ), by = svId])
    }
}
#=====================================================================================

#=====================================================================================
# create the composite bin coverage file
#-------------------------------------------------------------------------------------
message("printing results")
saveRDS(
    bins, 
    file = paste(env$COVERAGE_PREFIX, 'rds', sep = ".")
)
if(env$ASSESS_CNV_BINS) saveRDS(
    cnvBins, 
    file = paste(env$COVERAGE_PREFIX, 'cnvs.rds', sep = ".")
)
#=====================================================================================
