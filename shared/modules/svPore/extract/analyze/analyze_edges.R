#=====================================================================================
# definition of terms used in svPore analyze and find
#-------------------------------------------------------------------------------------
# read          a single, contiguous DNA molecule sequenced by a nanopore, as assessed by the basecaller (could be chimeric!)
# node          a specific numbered position on a genome strand; two nodes define an edge
# edge          the connection between two nodes, either an alignment or an SV junction
# alignment     an edge that corresponds to a portion of a read aligned contiguously to the genome by minimap2
# junction      an edge that nominates an SV by connecting two distant alignments
#                   each read has a series of edges as: alignment[-junction-alignment...]
# readPath      a sequence of nodes and junction properties that completely describes one read (or segment thereof)
#                   a readPath conforms to: node-alignmentEdge-[-node-junctionEdge-node-alignmentEdge...]-node
#                   where junctionEdge carries insertSize, i.e., information on the number of read/query bases at the junction
# flanks        the two alignments on either side of a junction
#                   reads are split into segments at junctions flanked by low quality alignments (among other reasons)
# block         a subset of a read aligned to a single strand of a single chromosome
#                   blocks may contain deletion, duplication and insertion junctions, but never inversions or translocations
# bandwidth     a base span on the reference genome strand minus the corresponding span on a read
#                   no SVs are called in a read span with bandwidth < MIN_SV_SIZE as they likely arise from misaligned low quality bases
#                   junctions that fail bandwidth do NOT split a read into segments
# segment       a subset of a readPath considered to be non-chimeric, arising from a single source DNA molecule
#                   unlike blocks, which are retained in a single segment row, segments are split from each other into separate rows
#                   a segment may contain SV junctions if they could be confirmed by multiple independent source molecules
# canonical     the strand orientation agreed to represent a given junction or readPath to allow comparison across reads and strands 
#                   junction canonicity is determined from the nodes flanking the junction
#                   read/segment canonicity is determined from the outermost alignment nodes
#-------------------------------------------------------------------------------------
# in total, in decreasing order of hierarchy:
#   all reads analyzed here contain at least two alignment edges and one junction edge (maybe more)
#   duplex reads, i.e., with the same readPath on opposite strands, are collapsed to single reads on just one strand
#   reads may be split to multiple segments at chimeric, low-quality, or unconfirmed junctions
#   a segment may contain one or more blocks separated by inversion or translocation junctions that may or may not call SVs
#   a block may contain one or more alignments separated by insertion, deletion, or duplication junctions that may or may not call SVs
#   an alignment may contain smaller indel operations in its CIGAR string that are never called as SVs
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)    
    library(e1071)  # provides the svm classifier
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
        'EXTRACT_STEP_DIR',
        'DATA_NAME',
        'EDGES_SV_FILE',
        'EDGES_TMP_FILE',
        'ANALYZE_PREFIX'
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
sourceScripts(env$EXTRACT_STEP_DIR, c(
    'edges','svm','duplex'
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
#=====================================================================================

# =====================================================================================
# initial loading and parsing of edges
# -------------------------------------------------------------------------------------
edges <- loadEdges("sv")
# edges <- edges[qName %in% c("714a3db0-c419-42be-b680-0efad4a9969a", "80b6332a-ca0b-4040-9549-b683edd73059")]
# edges[, cigar := "X"]
edges[, c(controlAdapterScores) := NULL]
setkey(edges, readI, blockN, edgeN)
edges <- checkReadBandwidth(edges)
#=====================================================================================

#=====================================================================================
# adapter splitting of chimeric molecules that derive from failure of the 
# basecaller to recognize that a new molecule had entered the pore
#-------------------------------------------------------------------------------------
if(is.null(env$SKIP_ADAPTER_CHECK)){
    trainingEdges <- loadEdges("tmp")
    trainingEdges[, cigar := NULL]
    trainingSet <- extractSvmTrainingSet(trainingEdges)
    svms <- trainAdapterClassifiers(trainingSet)
    rm(trainingEdges)
    adapterCheck <- checkJunctionsForAdapters(svms, edges[insertSize >= 5])
    edges <- updateEdgesForAdapters(edges, adapterCheck)
    rm(trainingSet, svms, adapterCheck)    
} else {
    message("skipping adapter check")
    edges <- updateEdgesForAdapters(edges)
}
#=====================================================================================

#=====================================================================================
# remove redundant duplex reads that were not fused during sequencing/basecalling
# ------------------------------------------------------------------------------------
if(is.null(env$SKIP_DUPLEX_CHECK)){
    reads <- collapseReads(edges, chromSizes)
    readMatches <- findMatchingReads(reads) 
    edges <- analyzeReadsNetwork(readMatches, reads, edges)
    rm(reads, readMatches)
} else {
    message("skipping duplex check")
}
#=====================================================================================

#=====================================================================================
# drop reads with no usable junctions and save for (multi-sample) SV finding
# ------------------------------------------------------------------------------------
edges <- dropReadsWithNoJunctions(edges, getMatchableJunctions)
saveRDS(
    edges, 
    paste(env$ANALYZE_PREFIX, "edges", "rds", sep = ".")
)
#=====================================================================================
