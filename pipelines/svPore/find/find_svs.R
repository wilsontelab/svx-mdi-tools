#=====================================================================================
# for definition of terms used in svPore find, see analyze/analyze_edges.R
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)    
    library(igraph) # provides graph/network analysis
    library(bit64)  # provides support for 64-bit integers
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'MODULES_DIR',
        'GENOMEX_MODULES_DIR',
        'ACTION_DIR',
        'GENOME',
        'GENOME_FASTA',
        'DATA_NAME',
        'EDGE_FILES',
        'FIND_PREFIX',
        'MANIFEST_PREFIX'
    ),
    integer = c(
        'N_CPU',
        'WINDOW_SIZE',
        'MIN_MAPQ',
        'MIN_ALIGNMENT_SIZE',
        'MIN_SV_SIZE',
        'JUNCTION_BANDWIDTH'
    ), 
    double = c(
        'MIN_ALIGNMENT_IDENTITY'
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general')) # , 'IUPAC', 'smith_waterman'
sourceScripts(file.path(rUtilDir, 'genome'),   c('chroms'))
sourceScripts(env$ACTION_DIR, c(
    'segments'
))
svPoreSharedDir <- file.path(env$MODULES_DIR, 'svPore')
sourceScripts(svPoreSharedDir, c(
    'constants','utilities','matching','filters','junctions'
))
setCanonicalChroms()
chromSizes <- loadChromSizes(env$WINDOW_SIZE)
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
# parse the project name and directory
if(env$FIND_MODE == "find"){
    projectName <- basename(env$OUTPUT_DIR)
    projectDir  <- env$OUTPUT_DIR
} else {
    projectName <- env$DATA_NAME
    projectDir  <- env$TASK_DIR
}
#=====================================================================================

# =====================================================================================
# initial loading and parsing of edges across potentially multiple samples/nanopores
# -------------------------------------------------------------------------------------
edges <- loadEdgesRds()
# edges[, cigar := "X"] ############### restrict tmp file size while developing
#=====================================================================================

#=====================================================================================
# matching individual SV junctions between molecules
# -------------------------------------------------------------------------------------
edges <- setMatchableFlag(edges)
edges <- setCanonicalNodes(edges, chromSizes)
junctionsToMatch <- getJunctionsToMatch(edges)
junctionHardCounts <- getJunctionHardCounts(junctionsToMatch)
junctionMatches <- findMatchingJunctions(junctionsToMatch) 
junctionsNetwork <- analyzeJunctionNetwork(junctionMatches, junctionHardCounts)
edges <- finalizeJunctionClustering(edges, junctionHardCounts, junctionsNetwork)
junctionClusters <- collapseJunctionClusters(edges)
rm(junctionsToMatch, junctionHardCounts, junctionMatches, junctionsNetwork)
#=====================================================================================

#=====================================================================================
# matching assembled segments to each other
# -------------------------------------------------------------------------------------
segmentEdges <- splitReadsToSegments(edges)
segments <- collapseSegments(segmentEdges, chromSizes)
nuclearSegments <- getNuclearSegments(segments)
segmentMatches <- findMatchingSegments(nuclearSegments)
segmentsNetwork <- analyzeSegmentsNetwork(segments, segmentMatches)
segmentClusters <- collapseSegmentClusters(segmentEdges, segmentsNetwork)
junctionClusters <- updateJunctionClustersForSegments(junctionClusters, segmentClusters)
#=====================================================================================

#=====================================================================================
# make a simple sample manifest and commit results to disk
# -------------------------------------------------------------------------------------
message("assembling windowCoverage files for data packages")    
coverages <- mergeWindowCoverageFiles()

message("creating sample manifest")    
manifest <- edges[, {
    jClusters <- junctionClusters[[sample]]
    isAlignment <- getAlignmentEdges(.SD)
    .(
        Project = projectName,
        Sample_ID = sample,
        Description = sample,
        Yield = as.numeric(coverages[[sample]]),
        nReads = length(unique(readI)),
        nEdges = .N,
        nJunctionEdges = sum(!isAlignment),
        nJunctionClusters = sum(jClusters > 0),
        nSingletonClusters = sum(jClusters == 1)
    )    
}, keyby = .(sample)]
outFile <- paste(env$MANIFEST_PREFIX, 'sample_manifest', 'csv', sep = ".")
fwrite(
    manifest, 
    file = outFile, 
    quote = FALSE, 
    sep = ",",    
    row.names = FALSE,   
    col.names = TRUE
)

# save chromosome metadata
message("saving chromosome metadata")   
chromosomeMetadata <- list(
    genome          = env$GENOME,
    canonicalChroms = canonicalChroms,
    chromIndex      = chromIndex,
    revChromIndex   = revChromIndex,
    chromSizes      = chromSizes
)
saveRDS(
    chromosomeMetadata, 
    paste(env$FIND_PREFIX, "chromosome_metadata", "rds", sep = ".")
)

# save results for the app
saveRDS(
    edges, 
    paste(env$FIND_PREFIX, "edges", "rds", sep = ".")
)
saveRDS(
    junctionClusters, 
    paste(env$FIND_PREFIX, "junctionClusters", "rds", sep = ".")
)
saveRDS(
    nuclearSegments, 
    paste(env$FIND_PREFIX, "segments", "rds", sep = ".")
)
saveRDS(
    segmentClusters, 
    paste(env$FIND_PREFIX, "segmentClusters", "rds", sep = ".")
)
#=====================================================================================
