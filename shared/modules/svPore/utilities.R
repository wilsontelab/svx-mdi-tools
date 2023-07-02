#-------------------------------------------------------------------------------------
# svPore analyze and find general support functions
#-------------------------------------------------------------------------------------

# load extract input data as edge tables
loadEdges <- function(type) {
    message(paste("loading", type, "edges"))
    x <- fread(
        switch(
            type,
            sv  = env$EDGES_SV_FILE,    # file with reads that had split alignments
            tmp = env$EDGES_TMP_FILE    # file with ~10K reads that did not have split alignments
        ),
        col.names  = edgesCols,
        colClasses = edgesColClasses,
        sep = "\t",
        quote = ""
    )
    x[, sample := env$DATA_NAME]
    x
}

# load and merge analyze edge tables across multiple samples
loadEdgesRds <- function(){
    EDGE_FILES <- strsplit(env$EDGE_FILES, "\\s+")[[1]]
    message(paste("loading sample edges", basename(EDGE_FILES[1])))
    edges <- readRDS(EDGE_FILES[1])
    if(length(EDGE_FILES) > 1){
        for(i in 2:length(EDGE_FILES)) {
            message(paste("loading sample edges", basename(EDGE_FILES[i])))
            edges <- rbind(edges, readRDS(EDGE_FILES[i]))
        }
    }
    setkey(edges, sample, readI, blockN, edgeN)
    edges
}

# merge window coverage files into a directory of bgz[.tbi] per sample
mergeWindowCoverageFiles <- function(){
    tmpDir <- file.path(env$TASK_DIR, "windowCoverage__WORKING")
    if(dir.exists(tmpDir)) unlink(tmpDir, recursive = TRUE, force = TRUE)
    dir.create(tmpDir, showWarnings = FALSE, recursive = FALSE)
    EDGE_FILES <- strsplit(env$EDGE_FILES, "\\s+")[[1]]
    coverages <- list()
    for(edgeFile in EDGE_FILES){
        coverageFile <- sub("analyze.edges.rds", "extract.windowCoverage.txt.bgz", edgeFile)
        indexFile <- paste(coverageFile, "tbi", sep = ".")
        file.copy(coverageFile, tmpDir)
        file.copy(indexFile, tmpDir)
        sample <- basename(dirname(edgeFile))
        coverage <- fread(
            cmd = paste("zcat", coverageFile),
            col.names  = c("chrom","pos","window","coverage"),
            colClasses = c("character","integer","integer","integer"),
            sep = "\t",
            quote = ""
        )
        coverages[[sample]] <- coverage[, median(coverage)]
    }
    zipFile <- paste(env$FIND_PREFIX, "windowCoverage", "zip", sep = ".")
    system2("zip", c("-jr", zipFile, tmpDir))
    unlink(tmpDir, recursive = TRUE, force = TRUE)
    coverages
}

# nodes are codified into an integer64 for streamlined comparison
# this function expands integer nodes out to chrom/strand/pos
parseSignedNodes <- function(chromSizes, nodes, side, canonical = FALSE) {
    genomeIs <- abs(nodes) # 1-referenced, per initialize_windows.pl
    chromIs <- Vectorize(function(i) which(chromSizes$nBasesThrough >= i)[1])(genomeIs) # sapply does not work with integer64!
    dt <- chromSizes[chromIs][, .(
        chrom = chrom, 
        chromIndex = chromIndex,
        refPos = as.integer(genomeIs - nBasesBefore),
        strand = ifelse(nodes > 0, "+", "-")
    )]    
    if(canonical) setnames(dt, c("cChrom","cChromIndex","cRefPos","cStrand"))
    setnames(dt, paste0(names(dt), side))
    dt
}
