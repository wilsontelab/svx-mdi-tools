# parse long-read molecules into validated and grouped SV paths

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
library(data.table)
library(e1071) # provides the svm classifier
library(parallel)
library(igraph)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOMEX_MODULES_DIR',
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
        'MIN_SV_SIZE',
        'MIN_MAPQ',
        'MIN_ALIGNMENT_SIZE'
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
sourceScripts(file.path(env$ACTION_DIR, 'analyze'), c('constants','utilities','svm','plot'))
setCanonicalChroms()
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

##############################
dir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svPore"
edgesRds            <- file.path(dir, "debug.edges.rds")
trainingSetRds      <- file.path(dir, "debug.trainingSet.rds")
svmsRds             <- file.path(dir, "debug.svms.rds")
jxnParametersRds    <- file.path(dir, "debug.jxnParameters.rds")
pathTypesTable      <- file.path(dir, "debug.pathTypes.txt")

#=====================================================================================
# initial loading and parsing of edges
#-------------------------------------------------------------------------------------
# edges <- loadEdges("sv")

# message("expanding edge metadata")
# edges <- cbind(edges, edges[, parseSignedWindow(node1, 1)], edges[, parseSignedWindow(node2, 2)])
# setkey(edges, qName)

# message("numbering edges and alignment blocks")
# edges[, ":="(
#     blockN = if(.N == 3){
#         if(edgeType[2] %in% splitTypes) 1:3 else 1
#     } else {
#         blocks <- rle(sapply(edgeType, "%in%", splitTypes))$lengths
#         unlist(sapply(1:length(blocks), function(i) rep(i, blocks[i]), simplify = FALSE))
#     },
#     edgeN = 1:.N, # number the edges in each molecule; alignments are odd, junctions are even
#     nEdges = .N,    
#     cigar = ifelse(is.na(cigar), cigar2, cigar),
#     passedFlankCheck = { # this DOES mask internal junctions, even if outermost molecule alignments are high quality
#         passed <- rep(FALSE, .N)
#         ai <- seq(1, .N, 2)
#         passed[ai] <- .SD[ai, 
#             mapQ >= env$MIN_MAPQ & 
#             eventSize >= env$MIN_ALIGNMENT_SIZE &
#             gapCompressedIdentity >= env$MIN_ALIGNMENT_IDENTITY
#         ]
#         ji <- seq(2, .N, 2)
#         passed[ji] <- sapply(ji, function(i) passed[i - 1] && passed[i + 1])
#         passed
#     }
# ), by = .(qName)]
# for(col in c("cigar2","blastIdentity2","gapCompressedIdentity2")) edges[[col]] <- NULL
# setkey(edges, qName, blockN, edgeN)

# ######### restrict tmp file size while developing
# edges[, cigar := NULL]

# # renderPlot("qualityDistribution", edges, suffix = "blastIdentity")
# # renderPlot("qualityDistribution", edges, suffix = "gapCompressedIdentity")
# # renderPlot("gapCompressedIdentity_vs_mapq", edges)

# message()
# str(edges)
# # print(edges[1:100])
# message("SAVING EDGES")
# saveRDS(edges, edgesRds)
# stop("XXXXXXXXXXXXXXXXXXXXX")
# message("reading RDS file")
# edges <- readRDS(edgesRds)
#=====================================================================================

#=====================================================================================
# calculate and apply additional edge quality metrics
#-------------------------------------------------------------------------------------
# message("checking indel bandwith to identify low quality alignment blocks")
# nodePairs <- list()
# fillNodePairs <- function(N){
#     i <- seq(2, N, 2)
#     x <- as.data.table(expand.grid(leftmost = i, rightmost = i))
#     x <- x[rightmost >= leftmost]
#     x[, nNodes := rightmost - leftmost + 1]
#     x[order(-nNodes)]
# } 
# checkBandwidth <- function(edges, i1, i2){
#     nQryBases <- edges[i2 + 1, xStart] - edges[i1 - 1, xEnd]
#     nRefBases <- edges[i2, xEnd] - edges[i1, xStart]
#     abs(nQryBases - nRefBases) >= env$MIN_SV_SIZE
# }
# edges[, passedBandwidth := {
#     if(.N == 1) TRUE # "A" blocks or "T|V" junctions 
#     else if(.N == 3) c(TRUE, checkBandwidth(.SD, 2, 2), TRUE) # single "D|U|I" junctions
#     else {
#         passed <- rep(TRUE, .N) # "ADAIADA" and other complex block paths
#         k <- as.character(.N)
#         if(is.null(nodePairs[[k]])) nodePairs[[k]] <<- fillNodePairs(.N)
#         for(j in 1:nrow(nodePairs[[k]])){
#             i1 <- nodePairs[[k]][j]$leftmost
#             i2 <- nodePairs[[k]][j]$rightmost
#             failed <- !checkBandwidth(.SD, i1, i2)
#             if(failed) passed[i1:i2] <- FALSE # don't write passes to prevent overwriting prior failures over wider junction spans
#         }
#         passed
#     }
# }, by = .(qName, blockN)]
# rm(nodePairs)

# message()
# str(edges)
# # print(edges[1:100])
# message("SAVING EDGES")
# saveRDS(edges, edgesRds)
# stop("XXXXXXXXXXXXXXXXXXXXX")
message("reading RDS file")
edges <- readRDS(edgesRds)

edges[, ":="(
    keptEdge = edgeType == edgeTypes$ALIGNMENT | (passedFlankCheck == TRUE & passedBandwidth == TRUE),
    keptJunction = edgeType != edgeTypes$ALIGNMENT & passedFlankCheck == TRUE & passedBandwidth == TRUE
)]
#=====================================================================================

#=====================================================================================
# perform adapter splitting
#-------------------------------------------------------------------------------------
# reads <- loadReads()
# trainingSet <- extractSvmTrainingSet()
# reads <- reads[molType == "J"]
# svms <- trainAdapterClassifiers(trainingSet)
# jxnParameters <- extractJunctionSvmParameters()
# rm(reads)

# saveRDS(trainingSet,    trainingSetRds)
# saveRDS(svms,           svmsRds)
# saveRDS(jxnParameters,  jxnParametersRds)

trainingSet     <- readRDS(trainingSetRds)
svms            <- readRDS(svmsRds)
jxnParameters   <- readRDS(jxnParametersRds)

adapterCheck <- checkJunctionsForAdapters(svms, jxnParameters)
edges <- updateEdgesForAdapters(edges, adapterCheck)
rm(trainingSet, svms, jxnParameters, adapterCheck)
#=====================================================================================

#=====================================================================================
message("scanning SV junctions for matching junctions in other molecules")
junctionsToMatch <- edges[keptJunction == TRUE, .SD, .SDcols = c("segmentName","blockN","edgeN","edgeType","chromIndex1","chromIndex2","node1","node2")]
junctionsToMatch[, chromPair := paste(sort(c(chromIndex1, chromIndex2)), collapse = ":"), by = c("segmentName","blockN","edgeN")]
edges <- merge(
    edges, 
    do.call(rbind, mclapply(1:nrow(junctionsToMatch), function(i){
        junctionsToMatch[i, getJunctionMatches(segmentName, blockN, edgeN, edgeType, chromPair, c(node1, node2))]
    }, mc.cores = env$N_CPU)), #
    by = c("segmentName","blockN","edgeN"), 
    all.x = TRUE
)
edges[, nMatchingSegments := sapply(matchingSegments, function(x) length(unlist(x)))]

message()
str(edges)
message()
str(junctionsToMatch)

# TODO: this had redundancy of IDENTICAL and DUPLEX
# need to process segments, then loop back to cleaned up edges
# ALSO: this doesn't obey the 1-window tolerance; must implement a shared name of matched nodes
gData <- edges[
    keptJunction == TRUE, 
    {
        isCanonical <- isCanonicalStrand(c(node1, node2))
        .(
            node1 = if(isCanonical) node1 else -node2,
            node2 = if(isCanonical) node2 else -node1,
            edgeType = edgeType
        )
    },
    by = c("segmentName","blockN","edgeN")
][, .(
    edgeTypes = paste(unique(edgeType), collapse = ","),
    nSegments = length(unique(segmentName))
), by = .(node1, node2)]

message()
str(gData)
message()
str(gData[grepl(",", edgeTypes)])
message()
str(gData[edgeTypes != "T" | nSegments > 1]) # a useful filter that discards interchromosomal singletons

# don't want to filter single edges for nSegments
# filter molecules groups for max(nSegments) or similar
# which means need to solve molecules groups

g <- graph_from_data_frame(gData[nSegments > 2], directed = TRUE)
message()
print(summary(g))
message()
str(components(g))

stop("XXXXXXXXXXXXXXXXXXXXXXX")

message("collapsing edges into source molecule segments, i.e., after adapter splitting")
segments <- edges[, .(
    nJxns = sum(edgeType != edgeTypes$ALIGNMENT),
    nKeptJxns = sum(keptJunction),
    pathType = paste0(ifelse(edgeType == edgeTypes$ALIGNMENT | keptJunction, edgeType, tolower(edgeType)), collapse = ""),
    isClosedPath = chromIndex1[1] == chromIndex2[.N] && strand1[1] == strand2[.N],
    nodePath = list(c(
        node1[1],
        c(rbind(node1[keptJunction], node2[keptJunction])),
        node2[.N]
    )),
    matchingSegmentsTmp = list(unique(unlist(matchingSegments)))
), by = .(segmentName)]
segments[, nMatchingSegments := sapply(matchingSegmentsTmp, function(x) length(unlist(x)))]

message()
str(segments)

message("scanning source molecule segments for matching junction patterns in other molecules")
segmentsToMatch <- segments[nMatchingSegments > 0]
segments <- merge(
    segments, 
    do.call(rbind, mclapply(1:nrow(segmentsToMatch), function(i){
        segmentsToMatch[i, getSegmentMatches(segmentName, unlist(matchingSegmentsTmp), unlist(nodePath))]
    }, mc.cores = env$N_CPU)),
    by = c("segmentName"), 
    all.x = TRUE
)
segments[, matchingSegmentsTmp := NULL]

# TODO: handle duplex and identical (collapse by removing replicates)
# TODO: assess closure of a set of overlapping segments

message()
str(segments)
message()
str(segments[nKeptJxns > 1])
# message()
# print(segments[, .N, by = .(nMatchingSegments)])
stop("XXXXXXXXXXXXXXXXXXXXX")


message()
str(edges)
# message("SAVING EDGES")
# saveRDS(edges, edgesRds)
# stop("XXXXXXXXXXXXXXXXXXXXX")
# message("reading RDS file")
# edges <- readRDS(edgesR
#=====================================================================================


stop("XXXXXXXXXXXXXXXXXXXXX")


qNames <- nodes[edgeType != edgeTypes$ALIGNMENT, unique(qName)]
saveRDS(nodes[qName %in% qNames], paste(env$DATA_FILE_PREFIX, "edges", "rds", sep = "."))












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

# message("collapsing nodes into segments")
# segments <- nodes[, .(
#     nJxns = sum(edgeType != edgeTypes$ALIGNMENT),
#     nKeptJxns = sum(edgeType != edgeTypes$ALIGNMENT & TRUE), # only count junctions that pass filters
#     junctions = NULL 
# ), by = .(segmentName)]


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