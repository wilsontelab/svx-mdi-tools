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
    library(yaml)
    library(stringr)
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
        'DATA_NAME',
        'BIN_SIZE'
    ), 
    integer = c(
        'N_CPU',
        'KMER_LENGTH',
        'N_ERRORS'
    )
))
env$ASSESS_JUNCTION_BINS <- Sys.getenv("ASSESS_JUNCTION_BINS") != "" # set by pipeline, not user
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
}#=====================================================================================

#=====================================================================================
# load the metadata across all samples
#-------------------------------------------------------------------------------------
message("loading sample metadata")
inFile <- paste(env$FIND_PREFIX, "metadata", "yml", sep=".")
metadata <- read_yaml(inFile)
metadata <- lapply(metadata, function(x) strsplit(as.character(x), "\\s+")[[1]])
#=====================================================================================

#=====================================================================================
# support functions for analyzing coverage
#-------------------------------------------------------------------------------------
getJunctionHalf <- function(sample, genRef, side){ # one base position error at junction point has ~no consequence on final value
    if(side == "L") substr(genRef, genRefIndexPos - jxnSideOffsets[[sample]], genRefIndexPos)
               else substr(genRef, genRefIndexPos, genRefIndexPos + jxnSideOffsets[[sample]])
}
getJunctionGC <- function(half1, half2){
    x <- paste0(half1, half2)
    (str_count(x, "G") + str_count(x, "C")) / str_length(x)
}
#-------------------------------------------------------------------------------------
getCompositeGenomes <- function(){
    binSizes <- as.integer(strsplit(env$BIN_SIZE, ",")[[1]])
    singleGenome <- data.table(
        genome  = env$GENOME,
        binSize = binSizes[1],
        tenBinsSize = binSizes[1] * 10,
        isComposite = FALSE
    )
    yamlFile <- file.path(env$GENOMES_DIR, env$GENOME, paste0(env$GENOME, ".yml"))
    if(!file.exists(yamlFile)) return(singleGenome)
    yaml <- read_yaml(yamlFile)
    if(is.null(yaml$composite) || !yaml$composite || is.null(yaml$compositeSources)) return(singleGenome)
    data.table(
        genome  = names(yaml$compositeSources),
        binSize = binSizes,
        tenBinsSize = binSizes * 10,
        isComposite = TRUE,
        compositeDelimiter = if(is.null(yaml$compositeDelimiter)) "_" else yaml$compositeDelimiter
    )
}
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
assignBins <- function(breakSpans, binSize){ # spans either on a whole chromosome, or adjacent to a specific junction
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
            start     = seq(bin1, bin2, binSize)
        ), by = .(breakStart)] else NULL
    )
    binAssignments[, end := start + binSize]
}
getSideAssignments <- function(genome, chrom_, breakSpans, 
                               jxn, lowPos, highPos, 
                               spanType = c("innerFlank", "outerFlank", "eventSpan")){
    bs <- breakSpans[breakStart %between% c(lowPos, highPos)]
    if(nrow(bs) == 0) return(NULL) 
    binAssignments <- assignBins(bs, genome$binSize) 
    binAssignments[, {
        overlaps <- pmin(breakEnd, end, highPos) - pmax(breakStart, start, lowPos)
        .(
            chrom = chrom_,
            coverage = weighted.mean(coverage, overlaps),
            svId = jxn$SV_ID,
            overlap = sum(overlaps),
            spanType = spanType
        )
    }, keyby = .(start, end)] 
}
#=====================================================================================

#=====================================================================================
# load and analyze junctions for copy number characterization
#-------------------------------------------------------------------------------------
if(env$ASSESS_JUNCTION_BINS){
    message("loading called and sequenced junctions")
    junctions <- readRDS(paste(env$FIND_PREFIX, 'structural_variants', 'rds', sep = "."))
    junctions <- junctions[N_SPLITS > 0] # we only analyze sequenced junctions with established breakpoint positions

    # tabulate the near-exact GC for the sequenced junction molecule over the "sweet spot" for reads (fast)
    # applies to all junctions and is independent of bins or genomes
    # later, intergenome translocations will typically be disregarded at GC bias profile will differ, with different net coverage, etc.
    message("calculating reference GC content for non-singleton junctions")
    genRefSize <- junctions[1, nchar(GEN_REF_1)]
    genRefHalfSize <- (genRefSize - 1) / 2
    genRefIndexPos <- genRefHalfSize + 1
    maxTLens <- as.integer(metadata$MAX_TLENS)
    jxnSideOffsets <- as.integer(maxTLens * 0.5) # the "sweet spot" where most read pairs would align
    names(jxnSideOffsets) <- metadata$SAMPLES
    jxnGC_null <- data.table(SV_ID = character(), sample = character(), gc = double())
    jxnGC_out <- do.call(rbind, lapply(metadata$SAMPLES, function(sample){ 
        sampleJxns <- junctions[junctions[[sample]] > 1] # ignore singletons, they aren't likely to be informativ re: CN
        if(nrow(sampleJxns) == 0) return(jxnGC_null)
        sampleJxns[, .(
            sample = sample,
            gc = getJunctionGC(
                getJunctionHalf(sample, GEN_REF_1, SIDE_1),
                getJunctionHalf(sample, GEN_REF_2, SIDE_2)
            )
        ), by = .(SV_ID)]   
    }))
}
#=====================================================================================

#=====================================================================================
# run bin coverage analysis per genome to allow:
#   1 - use of different bin sizes depending on genome size
#   2 - independent GC bias correction in downstream processes
#-------------------------------------------------------------------------------------
genomes <- getCompositeGenomes()
chunkSize <- 65536 # fixed by extract/base_coverage.pl
bins_out <- list()
# data.table(chrom = character(), start = integer(), end = integer(),
#                        gc = double(), excluded = integer(), genmap = double())
# for(sample in metadata$SAMPLES) bins_out[[sample]] <- double()
jxnBins_out <- data.table(svId = character(), spanType = character(), sample = character(), 
                          gc = character(), coverage = character(), overlap = character())
for(genomeI in 1:nrow(genomes)){
    genome <- genomes[genomeI]
    message(paste0(genome$genome, ", bin size = ", genome$binSize))
#=====================================================================================

#=====================================================================================
# load fixed-width bins data
#-------------------------------------------------------------------------------------
message("  loading bin metadata")
binsDir <- file.path(env$GENOMES_DIR, "bins", genome$genome, "fixed_width_bins")
if(!dir.exists(binsDir)) binsDir <- file.path(env$GENOMES_DIR, "../bins", genome$genome, "fixed_width_bins")
binsFile <- paste0(genome$genome, ".bins.size_", genome$binSize, ".k_", env$KMER_LENGTH, ".e_", env$N_ERRORS, ".bed.gz")
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
    message()
    message("WARNING: missing bins file")
    message(binsFile)
    message("proceeding with dummy values for gc, excluded, genmap")
    message("run the 'prepareBins' pipeline to create the missing file")
    message()
    setCanonicalChroms()
    chromSizes <- loadChromSizes(windowSize = genome$binSize)
    chromSizes[, .(
        start = (1:as.integer(nChromWindows) - 1) * genome$binSize,
        gc = 0.5,
        excluded = 0,
        genmap = 1 
    ), by = .(chrom)][, end := start + genome$binSize]
}
if(genome$isComposite) bins[, chrom := paste(chrom, genome$genome, sep = genome$compositeDelimiter)]
binChroms <- unique(bins$chrom)
#=====================================================================================

#=====================================================================================
# split junction candidates against the bins for CN estimation after GC bias correction in app
#-------------------------------------------------------------------------------------
message("  assembling allele and event coverage for non-singleton, same-chromosome junctions")

# initialize each sample
for(sample in metadata$SAMPLES){
    message(paste("    ", sample))
    samplePrefix  <- getSamplePrefix(sample)
    coverageIndex <- getSampleCoverageIndex(samplePrefix) 
    sampleJxns <- if(env$ASSESS_JUNCTION_BINS) junctions[
        JXN_TYPE != "T" &       # ignore translocations, flank and event spans aren't defined and genome bin sizes could differ
        junctions[[sample]] > 1 # ignore singletons for efficiency, they aren't likely to be informativ re: CN
    ] else NULL

    # initialize bin-level analysis of the junction and its flanks to establish allelic CN inside and outside of event
    offsets <- c(0, coverageIndex$cumNBreaks) 
    breaksFile <- paste(samplePrefix, env$GENOME, "coverage.breaks", sep = ".")
    countsFile <- paste(samplePrefix, env$GENOME, "coverage.counts", sep = ".")
    breaksFile <- file(breaksFile, "rb")
    countsFile <- file(countsFile, "rb")    
    binCoverage <- data.table(chrom = character(), start = integer(), end = integer(), coverage = double())
    jxnBinCoverage <- cbind(binCoverage, data.table(svId = character(), overlap = integer(), spanType = character()))
    nullJxnBinCoverage <- jxnBinCoverage

    # work on chromosome at a time
    coverageChroms <- unique(coverageIndex$chrom)
    coverageChroms <- coverageChroms[coverageChroms %in% binChroms]
    for(chrom_ in coverageChroms){ # this cannot be reliably parallelized, seemingly due to seek/read collisions?
        message(paste("      ", chrom_)) 

        # load the break spans, i.e., the contiguous segments of the genome covered by the same number of reads
        breakSpans <- coverageIndex[
            chrom == chrom_,
            {
                offset  <- offsets[.I]
                nBreaks <- cumNBreaks - offset
                seek(breaksFile, where = offset * 2, origin = "start", rw = "read")
                seek(countsFile, where = offset * 2, origin = "start", rw = "read")
                .(
                    breaks = readBin(breaksFile, "integer", n = nBreaks, size = 2, signed = FALSE) + start,
                    counts = readBin(countsFile, "integer", n = nBreaks, size = 2, signed = FALSE)
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
            bin1 = floor(breakStart / genome$binSize) * genome$binSize,
            bin2 = floor(breakEnd   / genome$binSize) * genome$binSize
        )]
        breakSpans[, ":="(sameBin = bin1 == bin2)]

        # calculate genome coverage by bin for this sample+chrom, used for coverage plots and CN segmentation
        binAssignments <- assignBins(breakSpans, genome$binSize)   
        binCoverage <<- rbind(binCoverage, binAssignments[, .(
            chrom = chrom_,
            coverage = weighted.mean(coverage, pmin(breakEnd, end) - pmax(breakStart, start))
        ), keyby = .(start, end)])

        # calculate GC per per event span four 10-bin junction flanks, for GC bias correction of junction CN in app
        if(!env$ASSESS_JUNCTION_BINS) next
        if(nrow(sampleJxns) == 0) next # again, already excludes singletons and unsequenced junctions
        setkey(breakSpans, breakStart)
        jxnBinCoverage <<- rbind(jxnBinCoverage, do.call(rbind, lapply(1:2, function(jxnSide){
            CHROM <- paste("CHROM", jxnSide, sep = "_")
            jxns <- sampleJxns[sampleJxns[[CHROM]] == chrom_]
            nJxns <- nrow(jxns)
            if(nJxns == 0) return(NULL)
            POS  <- paste("POS",  jxnSide, sep = "_")
            SIDE <- paste("SIDE", jxnSide, sep = "_")
            OTHER_POS <- paste("POS", if(jxnSide == 1) 2 else 1, sep = "_")
            do.call(rbind, mclapply(1:nJxns, function(i){  
                jxn <- jxns[i]
                lowPos <- min(jxn[[POS]], jxn[[OTHER_POS]])
                hghPos <- max(jxn[[POS]], jxn[[OTHER_POS]])
                rbind(
                    getSideAssignments(genome, chrom_, breakSpans, jxn, 
                                       lowPos - genome$tenBinsSize, lowPos, "outerFlank"), # the hosting allele CN
                    getSideAssignments(genome, chrom_, breakSpans, jxn, 
                                       hghPos, hghPos + genome$tenBinsSize, "outerFlank"),                 
                    getSideAssignments(genome, chrom_, breakSpans, jxn, 
                                       lowPos, min(lowPos + genome$tenBinsSize, hghPos), "innerFlank"), # the changing SV event DNA
                    getSideAssignments(genome, chrom_, breakSpans, jxn, 
                                       max(hghPos - genome$tenBinsSize, lowPos), hghPos, "innerFlank"),                 
                    if(jxn$SV_SIZE > 10e6) NULL # for now, skip very large, chromosome-scale events, slow and unlikely to be real
                    else getSideAssignments(genome, chrom_, breakSpans, jxn, 
                                            lowPos, hghPos, "eventSpan") # the entirety of the changing SV event DNA
                )
            }, mc.cores = env$N_CPU))
        })))
    }
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
    if(env$ASSESS_JUNCTION_BINS && nrow(sampleJxns) > 0){
        jxnBinCoverage <- merge(
            jxnBinCoverage, 
            bins[, .(chrom, start, gc, excluded)],
            by = c("chrom", "start"),
            all.x = TRUE
        )
        jxnBins_out <- rbind(jxnBins_out, jxnBinCoverage[excluded / genome$binSize < 0.1, .(
            sample   = sample,
            gc       = paste(round(gc, 3),       collapse = ","),
            coverage = paste(round(coverage, 1), collapse = ","),
            overlap  = paste(overlap,            collapse = ",")
        ), by = .(svId, spanType)])
    }
}
bins_out[[genome$genome]] <- bins
#=====================================================================================

#=====================================================================================
# close the genomes loop
#-------------------------------------------------------------------------------------
}
bins_out <- do.call(rbind, bins_out)
#=====================================================================================

#=====================================================================================
# create the composite bin coverage file
#-------------------------------------------------------------------------------------
message("printing results")
saveRDS(
    bins_out, 
    file = paste(env$COVERAGE_PREFIX, 'bins.rds', sep = ".")
)
if(env$ASSESS_JUNCTION_BINS) {
    saveRDS(
        jxnBins_out, 
        file = paste(env$COVERAGE_PREFIX, 'junctionBins.rds', sep = ".")
    )
    saveRDS(
        jxnGC_out, 
        file = paste(env$COVERAGE_PREFIX, 'junctionGC.rds', sep = ".")
    )
}
saveRDS(
    genomes, 
    file = paste(env$COVERAGE_PREFIX, 'genomes.rds', sep = ".")
)
#=====================================================================================
