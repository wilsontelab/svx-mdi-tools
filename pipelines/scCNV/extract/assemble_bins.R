#=====================================================================================
# merge bin coverage for all cells into a single 10x scCNV pipeline-compatible hdf5 file
# NB: this file is NOT compatible for use BY the 10x scCNV pipeline
#     but can be used INSTEAD of a true 10x scCNV output file by svx scCNV normalize
#-------------------------------------------------------------------------------------
#    https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/output/hdf5
#    also see accompanying summary file cell_ranger_file_format.txt
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
message("initializing")
suppressPackageStartupMessages({
    library(rhdf5)
    library(data.table)
    library(parallel)
    library(yaml)
})
#-------------------------------------------------------------------------------------
# parse environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
source(file.path(rUtilDir, 'utilities.R'))
checkEnvVars(list(
    string = c(
        'MDI_DIR',
        'GENOMEX_MODULES_DIR',
        'ALIGNMENT_DIR',
        'ALIGNMENT_FILES',
        'TASK_DIR',
        'DATA_NAME',
        'GENOME'
    ),
    integer = c(
        "N_CPU"
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'genome'), c('chroms'))
setCanonicalChroms()
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2)
#-------------------------------------------------------------------------------------
# set some constants
BIN_SIZE <- 20000L  # constant to match 10x scCNV pipeline
KMER_LENGTH <- 150L # TODO: expose as options?
N_ERRORS <- 1L
MIN_MAPQ <- 40L
#=====================================================================================

#=====================================================================================
# parse the genome-level information (genome_tracks)
#-------------------------------------------------------------------------------------
message("collecting genome-level information")
chroms <- unlist(revChromIndex)
chroms <- chroms[startsWith(chroms, "chr")]
#-------------------------------------------------------------------------------------
prepareBins <- paste( # as created by genomex-mdi-tools prepareBins
    env$GENOME, 
    'bins', 
    paste('size', BIN_SIZE,    sep = "_"), 
    paste('k',    KMER_LENGTH, sep = "_"), 
    paste('e',    N_ERRORS,    sep = "_"), 
    'bed',
    'gz',
    sep = "."
)
prepareBins <- file.path(env$MDI_DIR, 'resources', 'genomes', 'bins', env$GENOME, 'fixed_width_bins', prepareBins)
if(!file.exists(prepareBins)) stop(paste("file not found:", prepareBins))
prepareBins <- fread(prepareBins)
setnames(prepareBins, c(
    'chrom',
    'start',
    'end',
    'chrom_bin_i',
    'gc_fraction',
    'strand',
    'excluded',
    'gap',
    'bad',
    'umap',
    'genmap'
))
#-------------------------------------------------------------------------------------
pemap <- paste( # as created by genomex-mdi-tools prepareBins
    env$GENOME, 
    'pemap', 
    paste('r', KMER_LENGTH, sep = "_"), 
    paste('q', MIN_MAPQ,    sep = "_"), 
    paste('s', BIN_SIZE,    sep = "_"), 
    'bed',
    'gz',
    sep = "."
)
pemap <- file.path(env$MDI_DIR, 'resources', 'genomes', 'pemap', env$GENOME, pemap)
if(!file.exists(pemap)) stop(paste("file not found:", pemap))
pemap <- fread(pemap)
setnames(pemap, c(
    'chrom',
    'start',
    'end',
    'bin_i',
    'pemap',
    'strand'
))
prepareBins <- merge(
    prepareBins, 
    pemap[, .SD, .SDcols = c("chrom", "start", "end", "pemap")], 
    by = c("chrom", "start", "end"),
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
)
#-------------------------------------------------------------------------------------
prepareBins[, chrom_bin_i := floor(start / BIN_SIZE) + 1]
cols <- list(
    gc_fraction = 'gc_fraction',
    # mappability = 'genmap'
    mappability = 'pemap'
)
genome_tracks <- lapply(names(cols), function(colIn){
    colOut <- cols[[colIn]]
    x <- lapply(chroms, function(chrom_){
        chrom_bin_i <- prepareBins[chrom == chrom_, chrom_bin_i]
        valuesOut <- rep(0, max(chrom_bin_i))
        valuesOut[chrom_bin_i] <- unname(unlist(prepareBins[chrom == chrom_, ..colOut]))
        valuesOut
    })
    names(x) <- chroms
    x
})
names(genome_tracks) <- names(cols)
nChromBins <- unname(sapply(chroms, function(chrom_) prepareBins[chrom == chrom_, max(chrom_bin_i)]))
#=====================================================================================

#=====================================================================================
# parse the sample-level information (metadata)
#-------------------------------------------------------------------------------------
message("collecting sample-level information")
metadata <- list(
    assembly  = env$GENOME,
    sample_id = env$DATA_NAME,
    organism = switch(
        env$GENOME,
        hg38   = "Homo_sapiens",
        GRCh38 = "Homo_sapiens",
        mm10   = "Mus_musculus",
        GRCm38 = "Mus_musculus",
        NA
    ),
    pipeline = "svx-mdi-tools/scCNV"
)
#=====================================================================================

#=====================================================================================
# parse the cell-level information (per_cell_summary_metrics)
#-------------------------------------------------------------------------------------
message("collecting cell-level information")
per_cell_summary_metrics <- data.table(alignment_file = strsplit(env$ALIGNMENT_FILES, "\\s")[[1]])
per_cell_summary_metrics[, ':='(
    bins_file = paste(alignment_file, 'bins', 'gz', sep = "."),
    junctions_file = paste(alignment_file, 'junctions', 'gz', sep = "."),
    cell_id = .I - 1,    
    cell_name = sapply(alignment_file, function(x){
        x <- strsplit(basename(x), '\\.')[[1]]
        i <- length(x) - 1
        if(x[i] == "name") i <- i - 1
        x[i]
    }),    
    total_num_reads = sapply(alignment_file, function(x){ # all input read pairs from fastq
        read_yaml(paste(x, "counts", "yml", sep = "."))$nReadPairs
    }),    
    num_discarded_reads = sapply(alignment_file, function(x){           # input read pairs failing either mapping or mapq
        read_yaml(paste(x, "counts", "yml", sep = "."))$nDiscardedPairs # sum of 10x num_unmapped_reads + num_lowmapq_reads
    })
)]
constants <- list(
    bin_size = BIN_SIZE,
    chroms = chroms,
    genome_bins = as.integer(sum(nChromBins)),
    num_bins_per_chrom = as.integer(nChromBins),
    num_cells = as.integer(nrow(per_cell_summary_metrics)),
    num_chroms = as.integer(length(chroms))
)
#=====================================================================================

#=====================================================================================
# fill a table of bin counts over all cells
#-------------------------------------------------------------------------------------
message("creating bin-cell count matrices by chrom")
raw_counts <- lapply(per_cell_summary_metrics$bins_file, function(x){
    x <- fread(x)
    setnames(x, c('chromosome', 'chrom_bin_i', 'count'))    
    x[, chromosome := unlist(revChromIndex[chromosome])]
    x
})
names(raw_counts) <- per_cell_summary_metrics$cell_name
raw_counts <- lapply(chroms, function(chrom_){
    max_chrom_bin_i <- prepareBins[chrom == chrom_, max(chrom_bin_i)]
    nullCounts <- rep(0, max_chrom_bin_i)
    x <- sapply(per_cell_summary_metrics$cell_name, function(cellName){
        counts <- nullCounts
        rc <- raw_counts[[cellName]][raw_counts[[cellName]]$chromosome == chrom_]
        rc <- rc[chrom_bin_i <= max_chrom_bin_i]
        counts[rc$chrom_bin_i] <- rc$count
        counts
    })
    unname(x)
})
names(raw_counts) <- chroms
per_cell_summary_metrics$num_mapped_dedup_reads <- unname(sapply(1:constants$num_cells, function(j){
    sum( sapply(chroms, function(chrom) sum(raw_counts[[chrom]][, j])) ) # ~final number of deduplicated read pairs used in analysis
}))
per_cell_summary_metrics[, num_duplicate_reads := total_num_reads - num_discarded_reads - num_mapped_dedup_reads] 
#=====================================================================================

#=====================================================================================
# create the hdf5 output file
#-------------------------------------------------------------------------------------
message("writing final hdf5 file, cnv_data.h5")
h5FileName <- "cnv_data.h5"
h5File <- file.path(env$TASK_DIR, h5FileName)
unlink(h5File)
h5createFile(h5File)
for(group in c('constants', 'genome_tracks', 'metadata', 'per_cell_summary_metrics', 'raw_counts')){
    h5write(get(group), h5File, group)
}
#=====================================================================================
