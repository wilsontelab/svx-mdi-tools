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
        'PLOT_PREFIX',
        'PLOTS_DIR'
    ),
    integer = c(
        'N_CPU',
        'MIN_SV_SIZE'
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
nodesRds            <- file.path(dir, "debug.nodes.rds")
trainingSetRds      <- file.path(dir, "debug.trainingSet.rds")
svmsRds             <- file.path(dir, "debug.svms.rds")
jxnParametersRds    <- file.path(dir, "debug.jxnParameters.rds")

#=====================================================================================
message("  loading nodes")
nodes <- fread(
    env$NODES_FILE,
    col.names = nodesCols,
    colClasses = nodesColClasses,
    sep = "\t",
    quote = ""
)
nodes <- cbind(nodes, nodes[, parseSignedWindow(node1, 1)], nodes[, parseSignedWindow(node2, 2)])
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

# saveRDS(nodes,          nodesRds)
# nodes           <- readRDS(nodesRds)
# =====================================================================================

# =====================================================================================
message("  checking indel bandwith to identify false, low quality SV segments")
checkBW <- function(xStart, xEnd, eventSize, firstJ, lastJ){
    nQryBases <- xEnd[lastJ] - xStart[firstJ]
    nRefBases <- sum(eventSize) # [edgeType %in% incrementRefTypes]
    abs(nQryBases - nRefBases) >= env$MIN_SV_SIZE # TRUE if segment passes bandwith check   
}
checkBandwidth <- function(nodes, firstI, lastI, results){
    x <- if(firstI == lastI) TRUE else nodes[firstI:lastI, { # simple alignments pass, e.g., in ATATA
        passedOuter <- checkBW(xStart, xEnd, eventSize, 1, .N)
        if(.N == 3 || !passedOuter) rep(passedOuter, .N) else { # if a complex junction sequence passed on outermost aligns, check for false sub-segments
            alnJ <- seq(1, .N, 2)
            subSegments <- as.data.table(expand.grid(fJ = alnJ, lJ = alnJ))[lJ > fJ & !(fJ == 1 & lJ == .N)]
            failedSubSegmentsK <- which(!apply(subSegments, 1, function(J) checkBW(xStart, xEnd, eventSize, J[1], J[2])))
            p <- rep(TRUE, .N)
            for(k in failedSubSegmentsK) p[subSegments[k, fJ:lJ]] <- FALSE
            p
        }    
    }]
    c(results, x)
}
nodes[, passedBandwidth := {
    if(.N == 1) TRUE else {
        firstI <- 1
        lastI <- 0
        results <- logical()
        for(i in 1:.N){
            if(edgeType[i] == edgeTypes$ALIGNMENT){ # segments are chains of D|I junctions ending in flanking alignments
                lastI <- i
            } else if(!(edgeType[i] %in% inlineSvTypes)){ # skip over D and I operations
                results <- c(checkBandwidth(.SD, firstI, lastI, results), TRUE) # T, V, U operations are always passed
                firstI <- i + 1
            }
        }
        checkBandwidth(.SD, firstI, lastI, results)
    }
}, by = "qName"]

saveRDS(nodes,          nodesRds)
# nodes           <- readRDS(nodesRds)
#=====================================================================================

#=====================================================================================
message("  running adapter splitting using support vector machines")
message("    loading reads")
reads <- fread(
    env$SEQUENCES_FILE,
    col.names = readsCols,
    colClasses = readsColClasses,
    sep = "\t",
    quote = ""
)
trainingSet <- extractSvmTrainingSet(nodes, reads)
svms <- trainAdapterClassifiers(trainingSet)

saveRDS(trainingSet,    trainingSetRds)
saveRDS(svms,           svmsRds)
# trainingSet     <- readRDS(trainingSetRds)
# svms            <- readRDS(svmsRds)

jxnParameters <- extractJunctionSvmParameters(nodes, reads)
rm(reads)

saveRDS(jxnParameters,  jxnParametersRds)
# jxnParameters   <- readRDS(jxnParametersRds)

adapterCheck <- checkJunctionsForAdapters(svms, jxnParameters)
nodes <- updateNodesForAdapters(nodes, adapterCheck)

######################
# nodes <- splitChimericMolecules(nodes)

rm(adapterCheck)

# saveRDS(nodes,          nodesRds)
# nodes           <- readRDS(nodesRds)

qNames <- nodes[edgeType != edgeTypes$ALIGNMENT, unique(qName)]
saveRDS(nodes[qName %in% qNames], paste(env$DATA_FILE_PREFIX, "svNodes", "rds", sep = "."))

#=====================================================================================

#=====================================================================================
# message("  aggregating junction calls")
# # TODO: apply SV filters prior to aggregating junctions?
# junctions <- nodes[edgeType != edgeTypes$ALIGNMENT, .(
#     nInstances = .N, # how many times this junction was encountered in any path
#     fractionChimeric = sum(hasAdapter) / .N, # how many of those times it was chimeric, i.e., has at least one found adapter
#     fractionPassedBandwidth = sum(passedBandwidth) / .N,
#     edgeType = edgeType[1]

#     # additional quality metrics
#     # lists of molecule and segment ids

# ), by = .(junction)]

# message("  collapsing nodes into segments")
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