#=====================================================================================
# assemble the parts of an initial BioConductor SummarizedExperiment-like object set
# although we do NOT use a SE because we end up with different binning per cell
#-------------------------------------------------------------------------------------
message("parsing genome bin descriptions")
autosomes <- paste0("chr", 1:100)
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
        autosome    = chrom %in% autosomes,
        chrom_bin_n = 1:nBins  
    )
    badRegions <- badRegions[chrom == constants$chroms[i], 2:3]
    setkey(badRegions, start, end)     
    dt[, bad_region := !is.na(foverlaps(dt, badRegions, type = "any", mult = "first", which = TRUE))]
    dt
}, mc.cores = env$N_CPU))
rowRanges[, bin_n := 1:.N]
rm(autosomes, badRegions, genome_tracks)
# -------------------------------------------------------------------------------------
message("parsing cell metrics")
cell_ids <- as.character(per_cell_summary_metrics$cell_id)
colData <- as.data.table(per_cell_summary_metrics)
colData[, cell_id := cell_ids]
setkey(colData, "cell_id")
rm(per_cell_summary_metrics)
# -------------------------------------------------------------------------------------
message("parsing cell raw bin counts")
raw_counts <- do.call(rbind, mclapply(constants$chroms, function(chrom){
    as.data.table(raw_counts[[chrom]][, 1:constants$num_cells]) # rows = all bins, columns = all cells
}, mc.cores = env$N_CPU))
setnames(raw_counts, cell_ids)
#=====================================================================================

#=====================================================================================
# permanently remove bins in bad genome regions; windows will span them
#-------------------------------------------------------------------------------------
message("removing bins in bad genome regions")
bad_bin_n  <-  rowRanges[bad_region == TRUE, bin_n]
raw_counts <- raw_counts[rowRanges$bad_region == FALSE]
rowRanges  <-  rowRanges[bad_region == FALSE]
#=====================================================================================

#=====================================================================================
# pre-aggregate the bins for all required window sizes
#-------------------------------------------------------------------------------------
message("pre-aggregating bins for all required window sizes")
window_sizes <- seq(1, env$MAX_WINDOW_BINS, 2)
windows <- mclapply(window_sizes, function(window_size){
    rowRanges[, {
        chrom_window_id <- floor((1:.N - 1) / window_size) + 1
        .SD[, .(
            start = min(start),
            end   = max(end),
            gc_fraction = mean(gc_fraction),
            mappability = mean(mappability),
            autosome = autosome[1]
        ), by = chrom_window_id]
    }, by = "chrom"]
}, mc.cores = env$N_CPU)
names(windows) <- paste("w", window_sizes, sep = "_")
#=====================================================================================
