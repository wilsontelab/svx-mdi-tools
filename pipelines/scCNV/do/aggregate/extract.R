# extract and reformat data from 10x Cell Ranger DNA input
# https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/output/hdf5

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
message("initializing")
suppressPackageStartupMessages({
    library(rhdf5)
    library(data.table)
    # library(SummarizedExperiment)
    # library(GenomicRanges)
    library(zoo)
    library(parallel)
})
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'INPUT_DIR',
        'INPUT_NAME',
        'GENOMES_DIR',
        'GENOME',
        'BAD_REGIONS_FILE',
        'TASK_DIR',
        'OUTPUT_DIR',
        'DATA_NAME',
        'PLOTS_DIR',
        'PLOT_PREFIX'
    ),
    integer = c(
        "PLOIDY",
        "MAX_WINDOW_BINS",
        "N_CPU"
    ),
    double = c(
        "N_SD_HALFCN",
        "MIN_MAPPABILITY",
        "TRANSITION_PROBABILITY"
    )
))
dir.create(env$PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

rDataFile <- file.path(env$OUTPUT_DIR, "test.RData")

#-------------------------------------------------------------------------------------
# source R scripts
classDir <- file.path(env$MODULES_DIR, 'classes/R/nbinomCountsGC')
sourceScripts(classDir, c('nbinomCountsGC_class', 'nbinomCountsGC_methods'))
classDir <- file.path(env$MODULES_DIR, 'classes/R/hmmEPTable')
sourceScripts(classDir, c('hmmEPTable_class', 'hmmEPTable_methods'))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# set some constants
autosomes <- paste0("chr", 1:100)
ploidyFactor <- (env$PLOIDY + 0.5) / env$PLOIDY
minBinCount <- (env$N_SD_HALFCN / (ploidyFactor - 1)) ** 2
#=====================================================================================

#=====================================================================================
# extract required data elements from Cell Ranger hdf5 file
# discard some Cell Ranger outputs, e.g., normalization and cell clusters
#-------------------------------------------------------------------------------------
message("loading Cell Ranger data")

# find the Cell Ranger hdf5 file, can be nested in a subfolder of inputDir/inputName
inputDir <- file.path(env$INPUT_DIR, env$INPUT_NAME)
h5FileName <- "cnv_data.h5"
h5Files <- list.files(
    path = inputDir, 
    pattern = h5FileName, 
    full.names = TRUE, 
    recursive = TRUE
)
h5File <- h5Files[1]
if(h5File == "") stop(paste("file", h5FileName, "not found in directory", inputDir))

# extract the required data elements
h5 = H5Fopen(h5File)
constants <- h5read(h5, "constants", bit64conversion = "int")
genome_tracks <- h5read(h5, "genome_tracks", bit64conversion = "int")
genome_tracks$is_mappable <- NULL
metadata <- h5read(h5, "metadata", bit64conversion = "int")
per_cell_summary_metrics <- h5read(h5, "per_cell_summary_metrics", bit64conversion = "int")
raw_counts <- h5read(h5, "raw_counts", bit64conversion = "int")
H5Fclose(h5)
#=====================================================================================

# # # save.image(file = rDataFile)
# # # stop("XXXXXXXXXXXXXXXXXXXX")
# # # load(rDataFile)

#=====================================================================================
# assemble the parts of an initial SummarizedExperiment
#-------------------------------------------------------------------------------------
message("parsing genome bin descriptions")

badRegions <- fread(env$BAD_REGIONS_FILE)[, 1:3]
setnames(badRegions, c("chrom", "start", "end"))
rowRanges <- do.call(rbind, mclapply(seq_along(constants$chroms), function(i){
    chrom <- constants$chroms[i]
    nBins <- constants$num_bins_per_chrom[i]
    start <- seq(from = 1, by = constants$bin_size, length.out = nBins)
    dt <- data.table(
        chrom = rep(chrom, nBins),
        start = as.integer(start),
        end   = as.integer(start + constants$bin_size - 1),
        gc_fraction = as.double(genome_tracks$gc_fraction[[chrom]]), 
        mappability = as.double(genome_tracks$mappability[[chrom]]),
        autosome    = chrom %in% autosomes
    )
    badRegions <- badRegions[chrom == constants$chroms[i], 2:3]
    setkey(badRegions, start, end)     
    dt[, bad_region := !is.na(foverlaps(dt, badRegions, type = "any", mult = "first", which = TRUE))]
    dt
}, mc.cores = env$N_CPU))

# # save.image(file = rDataFile)
# # stop("XXXXXXXXXXXXXXXXXXXX")
# # load(rDataFile)

# -------------------------------------------------------------------------------------
message("parsing cell metrics")
cell_ids <- as.character(per_cell_summary_metrics$cell_id)
colData <- as.data.table(per_cell_summary_metrics)
colData[, cell_id := cell_ids]
colData[, modal_CN := env$PLOIDY]
setkey(colData, "cell_id")
# rm(per_cell_summary_metrics)
# -------------------------------------------------------------------------------------
message("parsing cell raw bin counts")
raw_counts <- do.call(rbind, mclapply(constants$chroms, function(chrom){
    as.data.table(raw_counts[[chrom]][, 1:constants$num_cells]) # rows = all bins, columns = all cells
}, mc.cores = env$N_CPU))
setnames(raw_counts, cell_ids)
#=====================================================================================

# # save.image(file = rDataFile)
# # stop("XXXXXXXXXXXXXXXXXXXX")
# # load(rDataFile)

#=====================================================================================
# permanently remove bins in bad genome regions; windows will span them
#-------------------------------------------------------------------------------------
message("removing bins in bad genome regions")
raw_counts <- raw_counts[rowRanges$bad_region == FALSE]
rowRanges <- rowRanges[rowRanges$bad_region == FALSE]
#=====================================================================================

# save.image(file = rDataFile)
# stop("XXXXXXXXXXXXXXXXXXXX")
# load(rDataFile)

#=====================================================================================
# pre-calculate the bin corrections for all required window sizes
#-------------------------------------------------------------------------------------
message("pre-calculating bin corrections for all required window sizes")
window_sizes <- seq(1, env$MAX_WINDOW_BINS, 2)
rollingRanges <- mclapply(window_sizes, function(window_size){
    rmean <- function(x) rollmean( x, window_size, na.pad = TRUE, align = "center")
    # rmin  <- function(x) -rollmax(-x, window_size, na.pad = TRUE, align = "center")
    # rmax  <- function(x)  rollmax( x, window_size, na.pad = TRUE, align = "center")
    window_flank <- (window_size - 1) / 2
    windowUnit <- c(rep(FALSE, window_flank), TRUE, rep(FALSE, window_flank))
    rowRanges[,
        .(
            # start = as.integer(rmin(start)),
            # end   = as.integer(rmax(end)),
            gc_fraction = as.double(rmean(gc_fraction)), 
            mappability = as.double(rmean(mappability)),
            reference_window = rep(windowUnit, ceiling(.N / window_size))[1:.N]
        ),
        by = "chrom"
    ]
}, mc.cores = env$N_CPU)
names(rollingRanges) <- paste("w", window_sizes, sep = "_")
#=====================================================================================





save.image(file = rDataFile)
# stop("XXXXXXXXXXXXXXXXXXXX")
# load(rDataFile)

# #######################
# classDir <- file.path(env$MODULES_DIR, 'classes/R/nbinomCountsGC')
# sourceScripts(classDir, c('nbinomCountsGC_class', 'nbinomCountsGC_methods'))
# classDir <- file.path(env$MODULES_DIR, 'classes/R/hmmEPTable')
# sourceScripts(classDir, c('hmmEPTable_class', 'hmmEPTable_methods'))

#=====================================================================================
# characterize the individual cells
#-------------------------------------------------------------------------------------
message('characterizing individual cells')
windowBins <- rowRanges$autosome & rowRanges$mappability >= env$MIN_MAPPABILITY
makeGcBiasPlot <- function(cell_id, gc, NR_map, ER_gc){
    png(
        filename = paste(env$PLOT_PREFIX, "gc_bias", cell_id, "png", sep = "."),
        width = 1.5, 
        height = 1.5, 
        units = "in", 
        pointsize = 6,
        bg = "white",  
        res = 300,
        type = "cairo"
    )
    par(mar= c(4.1, 4.1, 0.1, 0.1))
    plot(gc, NR_map, xlim = c(0.35, 0.55), pch = 16, cex = 0.25, col = rgb(0, 0, 0, 0.2))
    points(gc, ER_gc, pch = 16, cex = 0.25, col = "red") # the fit negative binomial
    dev.off()
}
fitCell <- function(cell_id, minBinCount, pass){ # called twice: 1) learn about cell 2) optimize window size based on overdispersion

    # set the per-cell window size as the number of bins needed to obtain a median raw count >=minBinCount
    NR_raw_b <- raw_counts[[cell_id]][windowBins]
    NR_med_b <- median(NR_raw_b, na.rm = TRUE)
    if(is.na(NR_med_b) || NR_med_b == 0) return(NA)
    window_size <- ceiling(minBinCount / NR_med_b)
    window_size <- window_size + !(window_size %% 2) # make odd to ensure window centering

    # permanently remove cells with insufficient coverage, as determined by MAX_WINDOW_BINS
    if(window_size > env$MAX_WINDOW_BINS) return(NA)

    # calculate window sums and correct for mappability
    NR_raw_b <- raw_counts[[cell_id]] 
    rollingRanges <- rollingRanges[[paste("w", window_size, sep = "_")]]
    mappability <- rollingRanges[, ifelse(mappability < env$MIN_MAPPABILITY, NA, mappability)]
    NR_map_w <- unlist(sapply(constants$chrom, function(chrom){
        chromIs  <- rowRanges$chrom == chrom
        NR_raw_w <- rollsum(NR_raw_b[chromIs], window_size, na.pad = TRUE, align = "center")
        NR_raw_w / mappability[chromIs]
    }))

    # process required bin and cell information
    reference_windows <- rollingRanges$reference_window
    gc_w   <- rollingRanges$gc_fraction
    gc_wr  <- gc_w[reference_windows]
    gc_wra <- gc_w[reference_windows & rowRanges$autosome]
    NR_map_wr  <- NR_map_w[reference_windows] # THESE CAN BE NA DUE TO MABBABILITY
    NR_map_wra <- NR_map_w[reference_windows & rowRanges$autosome]
    chroms_wr <- rollingRanges$chrom[reference_windows]
    autosomes <- rollingRanges$chrom[reference_windows & rowRanges$autosome]

    # perform an initial GC bias fit that assumes the same CN for all autosomes
    label <- NULL # paste("cell", cell_id, "pass", pass, "fit", 1)
    fit <- new_nbinomCountsGC(NR_map_wra, gc_wra, binCN = env$PLOIDY, label = label)

    # use the initial fit to solve an initial CN estimate for all autosomes
    hmm <- viterbi(fit, NR_map_wra, gc_wra, asRle = FALSE,
                   chroms = autosomes, transProb = env$TRANSITION_PROBABILITY)
    cn_estimate <- ifelse(hmm$cn == hmm$maxCN, NA, hmm$cn) # bins at maxCN are unreliable as they might be >maxCN

    # revise to a final GC bias fit using the initial copy number estimates
    label <- NULL # paste("cell", cell_id, "pass", pass, "fit", 2)
    fit <- new_nbinomCountsGC(NR_map_wra, gc_wra, binCN = cn_estimate, label = label)

    # solve a final CN estimate for all chromosomes (not just autosomes)
    if(pass == 1){
        hmm <- viterbi(fit, NR_map_wr, gc_wr, asRle = FALSE, 
                       chroms = chroms_wr, transProb = env$TRANSITION_PROBABILITY)
        ER_gc <- predict(fit, gc_wr, type = 'adjustedPeak') * env$PLOIDY # use peak for visualization, unless it is unreliable
        cn <- NR_map_wr / ER_gc * env$PLOIDY
        makeGcBiasPlot(cell_id, gc_wr, NR_map_wr, ER_gc)        
    } else {
        hmm <- viterbi(fit, NR_map_w, gc_w, asRle = FALSE, 
                       chroms = rollingRanges$chrom, transProb = env$TRANSITION_PROBABILITY)
        ER_gc <- predict(fit, gc_w, type = 'adjustedPeak') * env$PLOIDY # use peak for visualization, unless it is unreliable
        cn <- NR_map_w / ER_gc * env$PLOIDY
    }

    # return our results
    list(
        rejected = TRUE, # caller must set to FALSE if keeping the cell
        window_size = window_size,
        NR_map = NR_map_w, # all overlapping moving windows
        fit = fit, # based on non-overlapping autosomal reference windows
        hmm = hmm, # in final pass, includes same windows as NR_map 
        ER_gc = ER_gc,
        cn = cn
    )
}
cellFits <- mclapply(cell_ids, function(cell_id){
    
    # first pass fit at the sensitivity expected for Poisson without over/under-dispersion
    x1 <- fitCell(cell_id, minBinCount, pass = 1)
    if(!is.list(x1)) return(NA)

    # second pass fit at a sensitivity adjusted for the specific cell's dispersion
    dcn <- x1$cn[x1$hmm$cn == env$PLOIDY] - env$PLOIDY
    scaleFactor <- env$N_SD_HALFCN * sd(dcn, na.rm = TRUE) / 0.5
    x <- fitCell(cell_id, minBinCount * scaleFactor, pass = 2)
    if(!is.list(x)) return(x1)
    x$rejected <- FALSE
    x

    # plot(1:length(x$cn), x$cn, typ = "n", 
    #      xaxs = "i", ylim = c(0, 4), xlab= NULL, ylab = "CN")
    # abline(h = 1:4, col = "grey")
    # abline(v = x$chromEnds, col = "grey")
    # points(1:length(x$cn), x$cn, pch = ".")
    # points(1:length(x$hmm$cn), x$hmm$cn, pch = 16, cex = 0.5, col = "red")
}, mc.cores = env$N_CPU)
names(cellFits) <- cell_ids

# assemble and organize the per-cell data
colData[, ':='(
    rejected    = sapply(cellFits, function(x) 
        if(is.list(x)) x$rejected else TRUE),
    window_size = sapply(cellFits, function(x) 
        if(is.list(x) && !is.null(x$window_size)) x$window_size else NA)
)]

str(colData)

# save.image(file = rDataFile)
# # stop("XXXXXXXXXXXXXXXXXXXX")
# # load(rDataFile)

# # str(cellFits)
# x <- sapply(cellFits, function(x) x$fit$theta)
# names(x) <- cell_ids
# print(x)
stop("XXXXXXXXXXXXXXXXXXXX")


createBatchEffectList <- function(keyObject, mse, gcBias, sampleIds, ...){

    # identify the good samples that have propagated forward in a specific sampleSet
    allSampleIds     <- names(gcBias$MCN_s) # all samples in the source
    unfailedSamples  <- sapply(gcBias$fit, is.list)
    batchedSampleJs  <- which(gcBias$batchReject_s == CONSTANTS$batched & unfailedSamples)
    batchedSampleIds <- allSampleIds[batchedSampleJs]
    batchedSampleIds <- batchedSampleIds[batchedSampleIds %in% sampleIds]
    goodSampleJs     <- which(gcBias$batchReject_s < CONSTANTS$rejected & unfailedSamples)
    goodSampleIds    <- allSampleIds[goodSampleJs]
    goodSampleIds    <- goodSampleIds[goodSampleIds %in% sampleIds]
    
    # calculate the percentile of each window in each cell relative to its nbinom GC fit
    x <- list()
    gc_w <- getCol(mse, 'gc_b_ma')
    x$lowertail <- as.data.table( lapply(cell_ids, function(cell_id){ 
        cumprob(
            cellFits[[cell_id]]$fit,
            binCounts = cellFits[[cell_id]]$NR_map, # getCol(mse, sampleId, 'NR_map_b'),
            fractionGC = gc_w,
            binCN = cellFits[[cell_id]]$hmm$cn #inverse.rle( gcBias$hmm[[sampleId]] )
            #binCN = inverse.rle( gcBias$CN_b[[sampleId]] )
        ) 
    }) )
    setnames(x$lowertail, cell_ids)

    # calculate the median percentile for each from among the bins and samples declared to be trustworthy
    x$lowertail$median <- apply(x$lowertail[, .SD, .SDcols = batchedSampleIds], 1, median, na.rm = TRUE)
    
    # use the median percentile for each bin to establish the correction
    x$readPerAllele_peak <- as.data.table( lapply(cell_ids, function(cell_id){
        gcFit <- cellFits[[cell_id]]$fit
        rpa <- predict(gcFit, fractionGC = gc_w, type = 'adjustedPeak')
        theta <- gcFit$theta
        
        # this and method below are very nearly equivalent...
        #qnbinom(x$lowertail$median, size=theta, mu=rpa)
        
        excess <- qnbinom(x$lowertail$median, size = theta, mu = rpa) -
                  qnbinom(0.5,                size = theta, mu = rpa)        
        rpa + excess # thus, expect more reads for bins with consistently high nbinom percentile
    }) )
    setnames(x$readPerAllele_peak, goodSampleIds)
    
    # return the result
    x$sampleIds <- goodSampleIds
    x
}



save.image(file = rDataFile)
stop("XXXXXXXXXXXXXXXXXXXX")
load(rDataFile)

# fit an initial HMM based on the initial curve fit
message('xxxx')

# revise the negative binomial GC curve fit based with copy number optimization
message('xxxx')

# revise the HMM based on the update GC bias correction
message('xxxx')

png(
    filename = file.path(env$OUTPUT_DIR, "test.png"),
    width = 3, 
    height = 3, 
    units = "in", 
    pointsize = 9,
    bg = "white",  
    res = 600,
    type = "cairo"
)
plot(hist(colSums(raw_counts), breaks = 50))
dev.off()
stop("XXXXXXXXXXXXXXXXXXXX")

# print(colSums(raw_counts))
# print(per_cell_summary_metrics$num_mapped_dedup_reads)
# stop("XXXXXXXXXXXXXXXXXXXX")

#-------------------------------------------------------------------------------------
# rm(genome_tracks)
# rowRanges <- makeGRangesFromDataFrame(
#     rowRanges,
#     keep.extra.columns = TRUE,
#     ignore.strand = TRUE,
#     starts.in.df.are.0based = FALSE
# )
# names(rowRanges) <- as.character(1:constants$genome_bins)

row.names(colData) <- cell_ids

se <- SummarizedExperiment(
    rowRanges = rowRanges,
    colData = colData,    
    assays = list(raw_counts = raw_counts),
    metadata  = list(
        constants = constants,
        metadata = metadata
    )
)
str(se)
rm(
    constants, genome_tracks, metadata, per_cell_summary_metrics, 
    rowRanges, colData, raw_counts
)
#=====================================================================================

# 0  cell_barcodes            H5I_DATASET STRING 329
# 1  cnvs                     H5I_GROUP             
# 2  constants                H5I_GROUP             
# 3  genome_tracks            H5I_GROUP             
# 4  metadata                 H5I_GROUP             
# 5  normalization_parameters H5I_GROUP             
# 6  normalized_counts        H5I_GROUP             
# 7  per_cell_summary_metrics H5I_GROUP             
# 8  raw_counts               H5I_GROUP             
# 9  summary_metrics          H5I_GROUP             
# 10 tree                     H5I_GROUP 

# for tables based on the 10x clustering tree (e.g., raw_counts):
#   rows 0 to N-1 correspond to the single cell copy number calls
#   rows N to 2N-2 represent the groups of cells defined by the hierarchical clustering

# cell_barcodes
#  chr [1:329(1d)] "AAACGGGTCTCGTGGG-1" "AAAGATGAGATCTGCT-1" ...

# cnvs
# List of 21
#  $ chr1 : int [1:9774, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr10: int [1:6535, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr11: int [1:6105, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr12: int [1:6007, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr13: int [1:6022, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr14: int [1:6246, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr15: int [1:5203, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr16: int [1:4911, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr17: int [1:4750, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr18: int [1:4536, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr19: int [1:3072, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr2 : int [1:9106, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr3 : int [1:8002, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr4 : int [1:7826, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr5 : int [1:7592, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr6 : int [1:7487, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr7 : int [1:7273, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr8 : int [1:6471, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chr9 : int [1:6230, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chrX : int [1:8552, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...
#  $ chrY : int [1:4588, 1:657] -128 -128 -128 -128 -128 -128 -128 -128 -128 -128 ...

# constants
# List of 7
#  $ bin_size          :integer64 20000 
#  $ chroms            : chr [1:21(1d)] "chr1" "chr2" "chr3" "chr4" ...
#  $ genome_bins       :integer64 136288 
#  $ num_bins_per_chrom:integer64 [1:21] 9774 9106 8002 7826 7592 7487 7273 6471 ... 
#  $ num_cells         :integer64 329 
#  $ num_chroms        :integer64 21 
#  $ num_nodes         :integer64 657 

# genome_tracks
# List of 4
#  $ gc_fraction:List of 21
#   ..$ chr1 : num [1:9774(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr10: num [1:6535(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr11: num [1:6105(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr12: num [1:6007(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr13: num [1:6022(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr14: num [1:6246(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr15: num [1:5203(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr16: num [1:4911(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr17: num [1:4750(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr18: num [1:4536(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr19: num [1:3072(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr2 : num [1:9106(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr3 : num [1:8002(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr4 : num [1:7826(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr5 : num [1:7592(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr6 : num [1:7487(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr7 : num [1:7273(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr8 : num [1:6471(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr9 : num [1:6230(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chrX : num [1:8552(1d)] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chrY : num [1:4588(1d)] 0 0 0 0 0 ...
#  $ is_mappable:List of 21  ##### threshold 90%
#   ..$ chr1 : Factor[1:9774(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr10: Factor[1:6535(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr11: Factor[1:6105(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr12: Factor[1:6007(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr13: Factor[1:6022(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr14: Factor[1:6246(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr15: Factor[1:5203(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr16: Factor[1:4911(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr17: Factor[1:4750(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr18: Factor[1:4536(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr19: Factor[1:3072(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr2 : Factor[1:9106(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr3 : Factor[1:8002(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr4 : Factor[1:7826(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr5 : Factor[1:7592(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr6 : Factor[1:7487(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr7 : Factor[1:7273(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr8 : Factor[1:6471(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr9 : Factor[1:6230(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chrX : Factor[1:8552(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chrY : Factor[1:4588(1d)] w/ 2 levels "FALSE","TRUE": 1 1 1 1 1 1 1 1 1 1 ...
#  $ mappability:List of 21
#   ..$ chr1 : num [1:9774, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr10: num [1:6535, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr11: num [1:6105, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr12: num [1:6007, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr13: num [1:6022, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr14: num [1:6246, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr15: num [1:5203, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr16: num [1:4911, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr17: num [1:4750, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr18: num [1:4536, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr19: num [1:3072, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr2 : num [1:9106, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr3 : num [1:8002, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr4 : num [1:7826, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr5 : num [1:7592, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr6 : num [1:7487, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr7 : num [1:7273, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr8 : num [1:6471, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chr9 : num [1:6230, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chrX : num [1:8552, 1] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ chrY : num [1:4588, 1] 0 0 0 0 0 ...
#  $ n_fraction :List of 21
#   ..$ chr1 : num [1:9774(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr10: num [1:6535(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr11: num [1:6105(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr12: num [1:6007(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr13: num [1:6022(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr14: num [1:6246(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr15: num [1:5203(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr16: num [1:4911(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr17: num [1:4750(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr18: num [1:4536(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr19: num [1:3072(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr2 : num [1:9106(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr3 : num [1:8002(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr4 : num [1:7826(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr5 : num [1:7592(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr6 : num [1:7487(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr7 : num [1:7273(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr8 : num [1:6471(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chr9 : num [1:6230(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chrX : num [1:8552(1d)] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ chrY : num [1:4588(1d)] 1 1 1 1 1 0.5 0 0 0 0 ...

# metadata
# List of 9
#  $ annotation      : chr [1(1d)] "gencode.vM17"
#  $ assembly        : chr [1(1d)] "GRCm38"
#  $ library_id      :'data.frame':       1 obs. of  2 variables:
#   ..$ gem_group : int [1(1d)] 1
#   ..$ library_id: chr [1(1d)] "LibraryNotSpecified"
#  $ organism        : chr [1(1d)] "Mus_musculus"
#  $ pipeline        : chr [1(1d)] "cnv"
#  $ pipeline_version: chr [1(1d)] "1.1.0"
#  $ reference_path  : chr [1(1d)] "/nfs/turbo/agc-data/refs/OLD_STRUCTURE/archived_refs/cellranger-dna/refdata-GRCm38-1.0.0"
#  $ sample_desc     : chr [1(1d)] ""
#  $ sample_id       : chr [1(1d)] "Sample_6774-SA-1"

# normalization_parameters
# List of 3
#  $ linear   : num [1:657(1d)] -0.6 0 -9.04 0 -5.04 -0.4 0 -0.4 -0.6 0 ...
#  $ quadratic: num [1:657(1d)] 0 0 16.72 0 4.56 ...
#  $ scale    : num [1:657(1d)] 1 1 1 1 1 1 1 1 1 1 ...

# normalized_counts 

# per_cell_summary_metrics
# List of 19
#  $ barcode                    : chr [1:329(1d)] "AAACGGGTCTCGTGGG-1" "AAAGATGAGATCTGCT-1" "AAAGCAAAGGCAACAC-1" "AAAGCAATCTGTGGCG-1" ...
#  $ cell_id                    :integer64 [1:329] 0 1 2 3 4 5 6 7 ... 
#  $ effective_depth_of_coverage: num [1:329(1d)] 0.1249 0.0269 0.0945 0.0364 0.0769 ...
#  $ effective_reads_per_1Mbp   :integer64 [1:329] 659 125 455 160 397 692 659 668 ... 
#  $ est_cnv_resolution_mb      : num [1:329(1d)] 0.58 4.5 2.46 3.96 1.74 ...
#  $ frac_mapped_duplicates     : num [1:329(1d)] 0.139 0.126 0.136 0.121 0.135 ...
#  $ is_high_dimapd             :integer64 [1:329] 0 0 1 0 0 0 0 0 ... 
#  $ is_noisy                   :integer64 [1:329] 0 0 1 1 1 0 0 0 ... 
#  $ mean_ploidy                : num [1:329(1d)] 2.01 3.87 1.93 1.63 2.02 ...
#  $ normalized_dimapd          : num [1:329(1d)] 1.3 1.5 2.16 1.62 1.69 ...
#  $ normalized_mapd            : num [1:329(1d)] 0.0858 0.2547 0.1871 0.2436 0.1549 ...
#  $ num_duplicate_reads        :integer64 [1:329] 352062 68969 262703 91582 209077 361640 336197 351904 ... 
#  $ num_lowmapq_reads          :integer64 [1:329] 368158 136349 416266 225476 252830 379944 408515 390655 ... 
#  $ num_mapped_dedup_reads     :integer64 [1:329] 1797186 341058 1238990 437303 1081436 1887232 1794852 1821332 ... 
#  $ num_unmapped_reads         :integer64 [1:329] 8518 2210 8303 3357 5577 8324 10280 8205 ... 
#  $ ploidy_confidence          :integer64 [1:329] -2 -3 -4 -4 -4 -2 -2 -2 ... 
#  $ raw_dimapd                 : num [1:329(1d)] 1.3 1.5 2.16 1.62 1.69 ...
#  $ raw_mapd                   : num [1:329(1d)] 0.0858 0.2547 0.1871 0.2436 0.1549 ...
#  $ total_num_reads            :integer64 [1:329] 2525924 548586 1926262 757718 1548920 2637140 2549844 2572096 ... 

# raw_counts
# List of 21
#  $ chr1 : num [1:9774, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr10: num [1:6535, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr11: num [1:6105, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr12: num [1:6007, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr13: num [1:6022, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr14: num [1:6246, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr15: num [1:5203, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr16: num [1:4911, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr17: num [1:4750, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr18: num [1:4536, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr19: num [1:3072, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr2 : num [1:9106, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr3 : num [1:8002, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr4 : num [1:7826, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr5 : num [1:7592, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr6 : num [1:7487, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr7 : num [1:7273, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr8 : num [1:6471, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chr9 : num [1:6230, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chrX : num [1:8552, 1:657] 0 0 0 0 0 0 0 0 0 0 ...
#  $ chrY : num [1:4588, 1:657] 0 0 0 0 0 0 0 0 0 0 ...

# summary_metrics
# List of 2
#  $ 1  :List of 43
#   ..$ bad_bc_count            :integer64 54716513 
#   ..$ bc_q20_bases            :integer64 7011595037 
#   ..$ bc_q20_bases_fract      : num 0.993
#   ..$ bc_q30_bases            :integer64 6882451779 
#   ..$ bc_q30_bases_fract      : num 0.974
#   ..$ bc_tot_bases            :integer64 7063053802 
#   ..$ clipped_base_fract      : num 0.0767
#   ..$ conf_map_fract          : num 0.768
#   ..$ correct_bc_rate         : num 0.876
#   ..$ iqr_insert_size         :integer64 132 
#   ..$ mapped_bases            :integer64 114374060587 
#   ..$ median_insert_size      :integer64 215 
#   ..$ num_conf_mapped         :integer64 677936170 
#   ..$ num_far_chimeras        :integer64 4378768 
#   ..$ num_pos_chimeras        :integer64 289956795 
#   ..$ num_reads               :integer64 882904340 
#   ..$ num_single_mapped       :integer64 7659416 
#   ..$ num_unmapped            :integer64 19245892 
#   ..$ r1_contains_N           :integer64 982760 
#   ..$ r1_len                  :integer64 135 
#   ..$ r1_q20_bases            :integer64 57118627201 
#   ..$ r1_q20_bases_fract      : num 0.958
#   ..$ r1_q30_bases            :integer64 54594965350 
#   ..$ r1_q30_bases_fract      : num 0.916
#   ..$ r1_tot_bases            :integer64 59594833189 
#   ..$ r2_contains_N           :integer64 841674 
#   ..$ r2_len                  :integer64 151 
#   ..$ r2_q20_bases            :integer64 63756033888 
#   ..$ r2_q20_bases_fract      : num 0.956
#   ..$ r2_q30_bases            :integer64 60579776126 
#   ..$ r2_q30_bases_fract      : num 0.909
#   ..$ r2_tot_bases            :integer64 66658409790 
#   ..$ raw_bc_on_whitelist     :integer64 800914680 
#   ..$ reads_containing_N_fract: num 0.00207
#   ..$ si_q20_bases            :integer64 0 
#   ..$ si_q20_bases_fract      : num NaN
#   ..$ si_q30_bases            :integer64 0 
#   ..$ si_q30_bases_fract      : num NaN
#   ..$ si_tot_bases            :integer64 0 
#   ..$ single_mapped_fract     : num 0.00868
#   ..$ soft_clipped_bases      :integer64 9678675705 
#   ..$ total_bases             :integer64 126255320620 
#   ..$ unmapped_fract          : num 0.0218
#  $ all:List of 43
#   ..$ bad_bc_count            :integer64 54716513 
#   ..$ bc_q20_bases            :integer64 7011595037 
#   ..$ bc_q20_bases_fract      : num 0.993
#   ..$ bc_q30_bases            :integer64 6882451779 
#   ..$ bc_q30_bases_fract      : num 0.974
#   ..$ bc_tot_bases            :integer64 7063053802 
#   ..$ clipped_base_fract      : num 0.0767
#   ..$ conf_map_fract          : num 0.768
#   ..$ correct_bc_rate         : num 0.876
#   ..$ iqr_insert_size         :integer64 132 
#   ..$ mapped_bases            :integer64 114374060587 
#   ..$ median_insert_size      :integer64 215 
#   ..$ num_conf_mapped         :integer64 677936170 
#   ..$ num_far_chimeras        :integer64 4378768 
#   ..$ num_pos_chimeras        :integer64 289956795 
#   ..$ num_reads               :integer64 882904340 
#   ..$ num_single_mapped       :integer64 7659416 
#   ..$ num_unmapped            :integer64 19245892 
#   ..$ r1_contains_N           :integer64 982760 
#   ..$ r1_len                  :integer64 135 
#   ..$ r1_q20_bases            :integer64 57118627201 
#   ..$ r1_q20_bases_fract      : num 0.958
#   ..$ r1_q30_bases            :integer64 54594965350 
#   ..$ r1_q30_bases_fract      : num 0.916
#   ..$ r1_tot_bases            :integer64 59594833189 
#   ..$ r2_contains_N           :integer64 841674 
#   ..$ r2_len                  :integer64 151 
#   ..$ r2_q20_bases            :integer64 63756033888 
#   ..$ r2_q20_bases_fract      : num 0.956
#   ..$ r2_q30_bases            :integer64 60579776126 
#   ..$ r2_q30_bases_fract      : num 0.909
#   ..$ r2_tot_bases            :integer64 66658409790 
#   ..$ raw_bc_on_whitelist     :integer64 800914680 
#   ..$ reads_containing_N_fract: num 0.00207
#   ..$ si_q20_bases            :integer64 0 
#   ..$ si_q20_bases_fract      : num NaN
#   ..$ si_q30_bases            :integer64 0 
#   ..$ si_q30_bases_fract      : num NaN
#   ..$ si_tot_bases            :integer64 0 
#   ..$ single_mapped_fract     : num 0.00868
#   ..$ soft_clipped_bases      :integer64 9678675705 
#   ..$ total_bases             :integer64 126255320620 
#   ..$ unmapped_fract          : num 0.0218

# tree
# List of 3
#  $ Z               : num [1:4, 1:328] 117 272 29 2 133 156 29 2 104 330 ...
#  $ heterogeneity   :List of 21
#   ..$ chr1 : num [1:9774, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr10: num [1:6535, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr11: num [1:6105, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr12: num [1:6007, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr13: num [1:6022, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr14: num [1:6246, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr15: num [1:5203, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr16: num [1:4911, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr17: num [1:4750, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr18: num [1:4536, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr19: num [1:3072, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr2 : num [1:9106, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr3 : num [1:8002, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr4 : num [1:7826, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr5 : num [1:7592, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr6 : num [1:7487, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr7 : num [1:7273, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr8 : num [1:6471, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chr9 : num [1:6230, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chrX : num [1:8552, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#   ..$ chrY : num [1:4588, 1:328] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
#  $ is_cell_in_group: int [1:329, 1:328] 0 0 0 0 0 0 0 0 0 0 ...
