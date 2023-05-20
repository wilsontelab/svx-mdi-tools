#=====================================================================================
# definition of terms used in svPore analyze
#-------------------------------------------------------------------------------------
# read          a single, contiguous DNA molecule sequenced by a nanopore, as assessed by the basecaller (could be chimeric!)
# node          a specific numbered position or window on a genome strand; two nodes define an alignment or SV junction
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
# bandwidth     a block's base span on the reference genome strand minus the span on the read
#                   no SVs are called in a block with bandwidth < MIN_SV_SIZE as they likely arise from misaligned low quality bases
#                   junctions that failed bandwidth never split a read into segments
# segment       a subset of a readPath considered to be non-chimeric, arising from a single source DNA molecule
#                   unlike blocks, which are retained in a single segment row, segments are split from each other into separate rows
#                   a segment may contain SV junctions if they could be confirmed by multiple independent source molecules
# canonical     the strand orientation agreed to represent a given junction or readPath to allow comparison across reads and strands 
#                   junction canonicity is determined from the nodes flanking the junction
#                   segment  canonicity is determined from the outermost alignment nodes
#-------------------------------------------------------------------------------------
# in total, in decreasing order of hierarchy:
#   all reads analyzed here contain at least two alignment edges and one junction edge (maybe more)
#   duplex reads, i.e., with the same readPath on opposite strands, are collapsed to single reads on just one strand
#   reads may be split to multiple segments at chimeric, low-quality, or unconfirmed junctions
#   a segment may contain one or more blocks separated by inversion or translocation junctions that always call SVs
#   a block may contain one or more alignments separated by insertion, deletion, or duplication junctions that may or may not call SVs
#   an alignment may contain smaller indel operations in its CIGAR string that are never called as SVs
#=====================================================================================

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(e1071) # provides the svm classifier
    library(parallel)
    library(igraph) 
    library(bit64)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOMEX_MODULES_DIR',
        'GENOME_FASTA',
        'DATA_FILE_PREFIX',
        'TMP_DIR_WRK',
        'OUTPUT_DIR',
        'DATA_NAME',
        'EDGES_TMP_FILE',
        'EDGES_SV_FILE',
        'SEQUENCES_FILE',
        'COVERAGE_FILE',
        'PLOT_PREFIX',
        'PLOTS_DIR'
    ),
    integer = c(
        'N_CPU',
        'WINDOW_SIZE',
        'MIN_SV_SIZE',
        'MIN_MAPQ',
        'MIN_ALIGNMENT_SIZE',
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
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general', 'IUPAC', 'smith_waterman'))
sourceScripts(file.path(rUtilDir, 'genome'),   c('chroms'))
sourceScripts(file.path(env$ACTION_DIR, 'analyze'), c(
    'constants','utilities','matching','filters',
    'edges','svm','junctions','duplex','segments',
    'plot'
))
setCanonicalChroms()
chromSizes <- loadChromSizes(env$WINDOW_SIZE)
chromWindows <- expandChromWindows(chromSizes)
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

##############################
dir <- env$TASK_DIR
edges1Rds            <- file.path(dir, "debug.edges1.rds")
trainingSetRds       <- file.path(dir, "debug.trainingSet.rds")
svmsRds              <- file.path(dir, "debug.svms.rds")
jxnParametersRds     <- file.path(dir, "debug.jxnParameters.rds")
edges2Rds            <- file.path(dir, "debug.edges2.rds")
edges3Rds            <- file.path(dir, "debug.edges3.rds")
edges4Rds            <- file.path(dir, "debug.edges4.rds")

# =====================================================================================
# initial loading and parsing of edges and associated read and block-level quality metrics
# -------------------------------------------------------------------------------------
# edges <- loadEdges("sv")
# edges <- parseEdgeMetadata(edges, chromWindows)
# edges <- checkIndelBandwidth(edges)

# message("saving edges1 RDS file")
# saveRDS(edges, edges1Rds)
# message("reading edges1 RDS file")
# edges <- readRDS(edges1Rds)
#=====================================================================================

#=====================================================================================
# adapter splitting of chimeric molecules that derive from failure of the 
# basecaller to recognize that a new molecule had entered the pore
#-------------------------------------------------------------------------------------
# reads <- loadReads()
# trainingSet <- extractSvmTrainingSet()
# reads <- reads[molType == "J"]
# svms <- trainAdapterClassifiers(trainingSet)
# jxnParameters <- extractJunctionSvmParameters()
# rm(reads)

# message("saving SVM RDS files")
# saveRDS(trainingSet,    trainingSetRds)
# saveRDS(svms,           svmsRds)
# saveRDS(jxnParameters,  jxnParametersRds)
# message("reading SVM RDS files")
# trainingSet     <- readRDS(trainingSetRds)
# svms            <- readRDS(svmsRds)
# jxnParameters   <- readRDS(jxnParametersRds)

# adapterCheck <- checkJunctionsForAdapters(svms, jxnParameters)
# edges <- updateEdgesForAdapters(edges, adapterCheck)
# rm(trainingSet, svms, jxnParameters, adapterCheck)

# message("saving edges2 RDS file")
# saveRDS(edges, edges2Rds)
# message("reading edges2 RDS file")
# edges <- readRDS(edges2Rds)
#=====================================================================================

#=====================================================================================
# matching individual SV junctions between molecules
# -------------------------------------------------------------------------------------
# edges <- dropReadsWithNoJunctions(edges)
# edges <- setCanonicalNodes(edges)
# junctionsToMatch <- getJunctionsToMatch(edges)
# junctionHardCounts <- getJunctionHardCounts(junctionsToMatch)
# junctionMatches <- findMatchingJunctions(junctionsToMatch) 
# junctionClusters <- analyzeJunctionNetwork(junctionMatches, junctionHardCounts)
# edges <- finalizeJunctionClustering(edges, junctionHardCounts, junctionClusters)
# rm(junctionsToMatch, junctionHardCounts, junctionMatches, junctionClusters)

# message("saving edges3 RDS file")
# saveRDS(edges, edges3Rds)
# message("reading edges3 RDS file")
# edges <- readRDS(edges3Rds)
#=====================================================================================

#=====================================================================================
# find and remove redundant duplex reads that were not fused during sequencing/basecalling
# -------------------------------------------------------------------------------------
# reads <- collapseReads(edges, chromSizes)
# duplexStatus <- findDuplexReads(reads)
# edges <- adjustEdgesForDuplex(edges, duplexStatus)
# rm(reads, duplexStatus)
# message("saving edges4 RDS file")
# saveRDS(edges, edges4Rds)
message("reading edges4 RDS file")
edges <- readRDS(edges4Rds)
#=====================================================================================

#=====================================================================================
# matching assembled segments to each other
# -------------------------------------------------------------------------------------
confirmedEdges <- splitReadsToSegments(edges)
segments <- collapseSegments(confirmedEdges, chromSizes)
nuclearSegments <- getNuclearSegments(segments)
segmentMatches <- findMatchingSegments(nuclearSegments)
segmentClusters <- analyzeSegmentsNetwork(nuclearSegments, segmentMatches)

message()
str(segmentClusters)
message()
print(segmentClusters[, .(nSegments = .N), by = .(indexSegment)][, .(nClusters = .N), by = .(nSegments)][order(nSegments)])


saveRDS(
    confirmedEdges, 
    paste(env$DATA_FILE_PREFIX, "edges", "rds", sep = ".")
)
saveRDS(
    nuclearSegments, 
    paste(env$DATA_FILE_PREFIX, "segments", "rds", sep = ".")
)
saveRDS(
    segmentClusters, 
    paste(env$DATA_FILE_PREFIX, "clusters", "rds", sep = ".")
)

# stop("XXXXXXXXXXXXXXXX")


# segmentPairs <- getSegmentPairs(segmentClusters)

# segments[, matchingSegmentsTmp := NULL]

# TODO: assess closure of a set of overlapping segments
#=====================================================================================

#=====================================================================================
# exploring single-instance junctions
# -------------------------------------------------------------------------------------
# singeltonJunctions <- pullSingletonJunctions(edges)
# stop("XXXXXXXXXXXXXXXX")

#=====================================================================================






#=====================================================================================
# message("aggregating junction calls")
# # TODO: apply SV filters prior to aggregating junctions?
# junctions <- nodes[edgeType != edgeTypes$ALIGNMENT, .(
#     nInstances = .N, # how many times this junction was encountered in any path
#     fractionChimeric = sum(hasAdapter) / .N, # how many of those times it was chimeric, i.e., has at least one found adapter
#     fractionPassedBandwidth = sum(passedBandwidth) / .N,
#     edgeType = edgeType[1]

#     # additional quality metrics
#     # lists of molecule and segment ids

# ), by = .(junction)]


# # molecules <- nodes[, .(
# #     EDGES = paste(EDGE_TYPE, collapse = ""),       
# #     N_EDGES = .N,
# #     MIN_MAPQ = min(MAPQ),
# #     MAX_MAPQ = max(MAPQ),
# #     MIN_ALN_SIZE = min(EVENT_SIZE[EDGE_TYPE == "A"]),
# #     MAX_ALN_SIZE = max(EVENT_SIZE[EDGE_TYPE == "A"]),    
# #     MIN_SV_SIZE = if(.N == 1) 0 else min(EVENT_SIZE[EDGE_TYPE != "A"]),
# #     MAX_SV_SIZE = if(.N == 1) 0 else max(EVENT_SIZE[EDGE_TYPE != "A"]),
# #     N_STRANDS = N_STRANDS[1],
# #     NODES = list(c(NODE1, NODE2[.N])),
# #     PATH = getPathSignature(c(NODE1, NODE2[.N]))  ## IS THIS CORRECT?
 
# # ), by = "QNAME"]
#=====================================================================================


#=====================================================================================
# find duplex paths sequenced as two independent segments (presumably rare)
# TODO: do this at the node level by marking duplicates?
# message("collapsing duplicate segment paths")
# molecules <- molecules[, .(
#     QNAME = QNAME[which.max(MIN_MAPQ)[1]],
#     EDGES = EDGES[1],
#     N_EDGES = N_EDGES[1],
#     N_JXNS = floor(N_EDGES[1] / 2),    
#     MIN_MAPQ = max(MIN_MAPQ),
#     MAX_MAPQ = max(MAX_MAPQ),
#     MIN_ALN_SIZE = MIN_ALN_SIZE[1],
#     MAX_ALN_SIZE = MAX_ALN_SIZE[1],
#     MIN_SV_SIZE = MIN_SV_SIZE[1],
#     MAX_SV_SIZE = MAX_SV_SIZE[1],
#     N_STRANDS = sum(N_STRANDS),
#     NODES = NODES[1]
# ), keyby = "PATH"]

# characterize composite SV paths, count molecules

#=====================================================================================

#=====================================================================================

# highQualitySvNodes <- nodes[edgeType != edgeTypes$ALIGNMENT & mapQ >= 50]
# trueSvs <- highQualitySvNodes[
#     passedBandwidth &
#     # (edgeType == edgeTypes$TRANSLOCATION | eventSize - abs(insertSize) >= 1000) & 
#     chrom1 != "chrM" & chrom2 != "chrM"
# ]
# notChimericSvs <- trueSvs[hasAdapter == FALSE]
# largeSvs <- trueSvs[edgeType == edgeTypes$TRANSLOCATION | abs(eventSize) >= env$MIN_SV_SIZE]
# largeSingletons <- largeSvs[nInstances == 1]
# largeMultiples  <- largeSvs[nInstances >= 3]
# junctions <- trueSvs[, .(
#     nInstances = .N, # how many times this junction was encountered in any path
#     fractionChimeric = sum(hasAdapter) / .N, # how many of those times it was chimeric, i.e., has at least one found adapter
#     edgeType = edgeType[1]
# ), by = .(junction)]
# largeJunctions <- largeSvs[, .(
#     nInstances = .N, # how many times this junction was encountered in any path
#     fractionChimeric = sum(hasAdapter) / .N, # how many of those times it was chimeric, i.e., has at least one found adapter
#     edgeType = edgeType[1]
# ), by = .(junction)]

# renderPlot("insertSize_vs_eventSize", highQualitySvNodes, suffix = "All")
# renderPlot("insertSize_vs_eventSize", trueSvs, suffix = "True")
# renderPlot("insertSize_vs_eventSize", notChimericSvs, suffix = "Achimeric")
# renderPlot("eventSizeDistribution", highQualitySvNodes, suffix = "All")
# renderPlot("eventSizeDistribution", trueSvs, suffix = "True")
# renderPlot("eventSizeDistribution", notChimericSvs, suffix = "Achimeric")

# renderPlot("fractionChimeric_vs_eventSize", highQualitySvNodes, suffix = "All")
# renderPlot("fractionChimeric_vs_eventSize", trueSvs, suffix = "True")

# renderPlot("fractionChimeric_vs_insertSize", highQualitySvNodes, suffix = "All")
# renderPlot("fractionChimeric_vs_insertSize", trueSvs, suffix = "True")
# renderPlot("fractionChimeric_vs_insertSize", largeSingletons, suffix = "Large-Singleton")

# stop("XXXXXXXXXXXXXXX")

# renderPlot("fractionChimeric_vs_nInstances", junctions, suffix = "True")
# renderPlot("fractionChimeric_vs_nInstances", largeJunctions, suffix = "Large")

# renderPlot("adapterScore_vs_position", suffix = "end5",   trainingSet, highQualitySvNodes, 5)
# renderPlot("adapterScore_vs_position", suffix = "start3", trainingSet, highQualitySvNodes, 3)

# noInsert.notChimeric <- largeSingletons[insertSize < 5  & hasAdapter == FALSE]
# noInsert.chimeric    <- largeSingletons[insertSize < 5  & hasAdapter == TRUE]
# insert.chimeric      <- largeSingletons[insertSize >= 5 & hasAdapter == TRUE]
# multiple.notChimeric <- largeMultiples[hasAdapter == FALSE]
# multiple.chimeric    <- largeMultiples[hasAdapter == TRUE]

# print(nrow(nodes))
# print(nrow(highQualitySvNodes))
# print(nrow(largeMultiples))
# print(nrow(multiple.notChimeric))
# print(nrow(multiple.chimeric))

# renderPlot("adapterScore_vs_position", suffix = "end5.noInsert.notChimeric",   trainingSet, noInsert.notChimeric, 5)
# renderPlot("adapterScore_vs_position", suffix = "end5.noInsert.chimeric",      trainingSet, noInsert.chimeric,    5)
# renderPlot("adapterScore_vs_position", suffix = "end5.insert.chimeric",        trainingSet, insert.chimeric,      5)
# renderPlot("adapterScore_vs_position", suffix = "end5.multiple.notChimeric",   trainingSet, multiple.notChimeric, 5)
# renderPlot("adapterScore_vs_position", suffix = "end5.multiple.chimeric",      trainingSet, multiple.chimeric,    5)

# renderPlot("adapterScore_vs_position", suffix = "start3.noInsert.notChimeric", trainingSet, noInsert.notChimeric, 3)
# renderPlot("adapterScore_vs_position", suffix = "start3.noInsert.chimeric",    trainingSet, noInsert.chimeric,    3)
# renderPlot("adapterScore_vs_position", suffix = "start3.insert.chimeric",      trainingSet, insert.chimeric,      3)
# renderPlot("adapterScore_vs_position", suffix = "start3.multiple.notChimeric",   trainingSet, multiple.notChimeric, 3)
# renderPlot("adapterScore_vs_position", suffix = "start3.multiple.chimeric",      trainingSet, multiple.chimeric,    3)

# rm(reads, trainingSet, jxnParameters, adapterCheck)
#=====================================================================================

# #=====================================================================================
# # static coverage plots by chromosome
# #-------------------------------------------------------------------------------------
# message("plotting coverage histogram")
# map <- fread(
#     cmd = paste("zcat", env$COVERAGE_FILE),
#     col.names = c(
#         "CHROM", # extract nodes fields
#         "WINDOW",
#         "COVERAGE"
#     ),
#     colClasses = c(
#         "character", # PAF fields
#         "integer",
#         "integer"
#     )
# )
# chroms <- map[, unique(CHROM)]
# nChroms <- length(chroms)
# chromIs <- as.list(1:nChroms)
# names(chromIs) <- chroms
# maxWindow <- map[, max(WINDOW)]
# peakCoverage <- map[COVERAGE > 0, getmode(COVERAGE)]
# ploidy <- 2
# readsPerAllele <- peakCoverage / ploidy
# maxY <- readsPerAllele * 5
# # png(
# #     # filename = paste(env$COVERAGE_FILE, "png", sep = "."),
# #     filename = paste("/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svLong/HISTOGRAM", "png", sep = "."),
# #     width = 3, 
# #     height = 4, 
# #     units = "in", 
# #     pointsize = 8,
# #     res = 600, 
# #     type = "cairo"
# # )
# # dt <- map[, .(COUNT = .N), by = COVERAGE][order(COVERAGE)]
# # plot(dt$COVERAGE, dt$COUNT, typ = "l", xlim = c(0, maxY))
# # abline(v = readsPerAllele * 0:maxY)
# # dev.off()
# #=====================================================================================
# message("plotting genome")
# png(
#     # filename = paste(env$COVERAGE_FILE, "png", sep = "."),
#     filename = paste("/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svPore/GENOME", "png", sep = "."),
#     width = 12, 
#     height = 0.5 * nChroms, 
#     units = "in", 
#     pointsize = 8,
#     res = 300, 
#     type = "cairo"
# )
# layout(matrix(1:nChroms, ncol = 1))
# map[, {
#     dt <- .SD[, .(
#         x = mean(WINDOW),
#         y = mean(COVERAGE)
#     ), by = floor(WINDOW / 1e4)]
#     par(mar = c(0.1, 4.1, 0.1, 0.1))
#     plot(dt$x, dt$y, typ="p", pch = 19, cex = 0.25, col = rgb(0, 0, 1, 0.25),
#          xaxt = "n",
#          xlab = paste(CHROM, "Coordinate (bp)"), ylab = "Coverage", 
#          xlim = c(1, maxWindow), ylim = c(0, maxY),
#          yaxs = "i")
#     abline(h = readsPerAllele * 0:maxY)
#     j <- goodJunctions[chromIndex1 == chromIs[[CHROM]]]
#     abline(v = j$windowIndex1 * 100, col = j$color)
#     j <- goodJunctions[chromIndex2 == chromIs[[CHROM]]]
#     abline(v = j$windowIndex2 * 100, col = j$color)
# }, by = CHROM]
# dev.off()