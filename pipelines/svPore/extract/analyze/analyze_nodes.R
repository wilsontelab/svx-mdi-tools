# parse long-read molecules into validated and grouped SV paths

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("  initializing")
library(data.table)
library(e1071) # provides the svm classifier
library(parallel)
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
        'NODES_FILE',
        'SEQUENCES_FILE',
        'COVERAGE_FILE',
        'PLOT_PREFIX'
    ),
    integer = c(
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general', 'IUPAC', 'smith_waterman'))
sourceScripts(file.path(env$ACTION_DIR, 'analyze'), c(
    'svm'
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# support functions
#-------------------------------------------------------------------------------------
isCanonicalStrand <- function(nodePair){
    o <- order(abs(nodePair))
    nodePair[o[1]] > 0
}
getPathSignature <- function(NODES){
    isCanonical <- isCanonicalStrand(c(NODES[1], NODES[length(NODES)]))
    if(!isCanonical) NODES <- -rev(NODES)
    paste(NODES, collapse = ":")
}
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
parseSignedWindow <- function(window) {
    strand <- ifelse(window > 0, "+", "-")
    window <- abs(window)
    chromIndex <- bitwShiftR(window, 24)
    data.table(
        chromIndex  = chromIndex,
        # chrom       = chromIndex ? $revChromIndex{$chromIndex} : "?",
        windowIndex = bitwAnd(window, 2**24 - 1),
        strand      = strand
    )
}
edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    MERGE_FAILURE = "M",
    PROPER        = "P",
    FUSED_MERGE_FAILURE = "F",
    REJECTED_INDEL = "R",
    FUSED_MERGE_FAILURE_REJECTED_INDEL = "Q"
)
#=====================================================================================

#=====================================================================================
message("  loading nodes")
nodes <- fread(
    env$NODES_FILE,
    col.names = c(
        "qName", # extract nodes fields
        "node1",
        "node2",
        "edgeType",
        "mapQ",
        "eventSize",
        "insertSize",
        "qStart",
        "qEnd",
        "nStrands"
    ),
    colClasses = c(
        "character", # PAF fields
        "integer",
        "integer",
        "character",
        "integer",
        "integer",
        "integer",
        "integer",
        "integer",
        "integer"
    )
)
message("  loading reads")
reads <- fread(
    env$SEQUENCES_FILE,
    col.names = c(
        "molType", # extract nodes fields
        "qName",
        "qSeq"
    ),
    colClasses = c(
        "character", # PAF fields
        "character",
        "character"
    )
)
#=====================================================================================

#=====================================================================================
message("  grouping and counting candidate SV junctions")
nodes[, ":="(
    edge = 1:.N # number the edges in each molecule; alignments are odd, junctions are even
), by = .(qName)]
nodes[edgeType != edgeTypes$ALIGNMENT, junction := {
    getPathSignature(c(node1, node2)) # oriented to the canonical strand for comparing molecules
}, by = .(qName, edge)]
nodes <- merge(
    nodes, 
    nodes[edgeType != edgeTypes$ALIGNMENT, .(
        nInstances = .N, # how many times this junction was encountered in any path
        nMolecules = length(unique(qName)) # how many unique molecules bore this junction (some may have crossed it more than once)
    ), by = .(junction)], 
    by = "junction", 
    all.x = TRUE
)
setkey(nodes, qName, edge)
#=====================================================================================

#=====================================================================================
message("  running adapter splitting using support vector machines")
trainingSet <- extractSvmTrainingSet(nodes, reads)
svms <- trainAdapterClassifiers(trainingSet)
jxnParameters <- extractJunctionSvmParameters(nodes, reads)
adapterCheck <- checkJunctionsForAdapters(svms, jxnParameters)
nodes <- merge(nodes, adapterCheck, by = c("qName","edge"), all.x = TRUE)
nodes[, isChimeric := hasAdapter3 | hasAdapter5]

dir <- "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svPore"
nodesRds <- file.path(dir, "debug.nodes.rds")
trainingRds <- file.path(dir, "debug.training.rds")
saveRDS(nodes, nodesRds)
saveRDS(trainingSet, trainingRds)
nodes <- readRDS(nodesRds)
trainingSet <- readRDS(trainingRds)

plotAdapterMatches(trainingSet, nodes, 5, "end")
plotAdapterMatches(trainingSet, nodes, 3, "start")

rm(reads, trainingSet, jxnParameters, adapterCheck)


# print(nodes[edgeType != edgeTypes$ALIGNMENT & mapQ == 60, .N, by = .(isChimeric)][order(isChimeric)])
# print(nodes[edgeType != edgeTypes$ALIGNMENT & mapQ == 60 & nInstances == 1, .N, by = .(isChimeric)][order(isChimeric)])
# print(nodes[edgeType != edgeTypes$ALIGNMENT & mapQ == 60 & nInstances == 2, .N, by = .(isChimeric)][order(isChimeric)])
# print(nodes[edgeType != edgeTypes$ALIGNMENT & mapQ == 60 & nInstances == 3, .N, by = .(isChimeric)][order(isChimeric)])
# print(nodes[edgeType != edgeTypes$ALIGNMENT & mapQ == 60 & nInstances == 4, .N, by = .(isChimeric)][order(isChimeric)])


# nodes[, , by = .(junction)]

# junctions <- nodes[edgeType != edgeTypes$ALIGNMENT & mapQ >= 50, .(
#     nInstances = .N, # how many times this junction was encountered in any path
#     fractionChimeric = sum(isChimeric) / .N # how many of those times it was chimeric, i.e., has at least one found adapter
# ), by = .(junction)]
# junctionCounts <- junctions[nInstances >= 3 & nInstances <= 15, .(
#     nJunctions = .N,
#     nInstances = mean(nInstances)
# ), by = .(fractionChimeric)]
# print(junctionCounts[order(fractionChimeric)])
# print(junctions[fractionChimeric == 0, .N, by = .(nInstances)][order(nInstances)])

# junctionCounts <- junctions[, .(
#     nJunctions = .N,
#     fractionChimeric = mean(fractionChimeric)
# ), by = .(nInstances)]
# # print(junctionCounts[order(nInstances)])


# pngFile <- file.path(dir, "debug.test.png")
# png(pngFile, width = 3, height = 3, units = "in", pointsize = 7, res = 600, type = "cairo")
# plot(
#     jitter(junctions$nInstances, amount = 0.5), 
#     jitter(junctions$fractionChimeric, amount = 0.1), 
#     xlab = "nInstances",
#     ylab = "fractionChimeric",
#     xlim = c(0, 30),
#     ylim = c(-0.1, 1.1),
#     pch = 19, cex = 0.25, 
#     col = rgb(0, 0, 1, 0.1)
# )
# dev.off()


# nodes <- merge(
#     nodes, 
#     nodes[edgeType != edgeTypes$ALIGNMENT, .(
#         fractionChimeric = sum(isChimeric) / .N # how frequently among nInstances this junction was chimeric, i.e., had at least one found adapter
#     ), by = .(junction)], 
#     by = "junction", 
#     all.x = TRUE
# )
# setkey(nodes, qName, edge)

# nodes[, segment := {
#     if(.N == 1 || !any(na.omit(isChimeric))) 1
#     else if(.N == 3) c(1, NA, 2)
#     else {
#         x <- 1
#         for(i in seq(2, .N, 2)){
#             x <- c(x, 
#                 if(isChimeric[i]) c(NA, x[i - 1] + 1)
#                 else rep(x[i - 1], 2)
#             )
#         }
#         x
#     }
# }, by = .(qName)]

# str(nodes)
# str(nodes[edgeType != "A"])
# str(nodes[!is.na(isChimeric)])
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






# message("counting junctions")
# junctions <- nodes[EDGE_TYPE != "A", .(
#     EDGE_TYPE = EDGE_TYPE[1],
#     N = .N,
#     MAPQ = max(MAPQ),
#     EVENT_SIZE = EVENT_SIZE[1]
# ), by = .(NODE1, NODE2)]
# goodJunctions <- junctions[MAPQ > 55 & N > 3 & (EVENT_SIZE > 1e4 | EDGE_TYPE == "T")]
# expanded1 <- parseSignedWindow(goodJunctions$NODE1)
# setnames(expanded1, c("chromIndex1", "windowIndex1", "strand1"))
# expanded2 <- parseSignedWindow(goodJunctions$NODE2)
# setnames(expanded2, c("chromIndex2", "windowIndex2", "strand2"))
# goodJunctions <- cbind(goodJunctions, expanded1, expanded2)
# goodJunctions[, color := ifelse(EDGE_TYPE == "T", "red3", "orange2")]


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
# #=====================================================================================

# #=====================================================================================
# # message("collapsing nodes into molecules")
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

# # # find duplex paths sequenced as two molecules (presumably rare)
# # message("collapsing duplicates")
# # molecules <- molecules[, .(
# #     QNAME = QNAME[which.max(MIN_MAPQ)[1]],
# #     EDGES = EDGES[1],
# #     N_EDGES = N_EDGES[1],
# #     N_JXNS = floor(N_EDGES[1] / 2),    
# #     MIN_MAPQ = max(MIN_MAPQ),
# #     MAX_MAPQ = max(MAX_MAPQ),
# #     MIN_ALN_SIZE = MIN_ALN_SIZE[1],
# #     MAX_ALN_SIZE = MAX_ALN_SIZE[1],
# #     MIN_SV_SIZE = MIN_SV_SIZE[1],
# #     MAX_SV_SIZE = MAX_SV_SIZE[1],
# #     N_STRANDS = sum(N_STRANDS),
# #     NODES = NODES[1]
# # ), keyby = "PATH"]

# # # characterize individual SV junctions, count molecules

# # # characterize composite SV paths, counts molecules

# # # above required strand-aware comparison, similar to NODE_PATH1 == reverse(-NODE_PATH2)

# #=====================================================================================
