# sort and group molecules based on the SVs they carry

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'OUTPUT_DIR',
        'DATA_NAME',
        'NODES_FILE',
        'COVERAGE_FILE'
    ),
    integer = c(
        'N_CPU'
    )
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
#=====================================================================================

#=====================================================================================
# static coverage plots by chromosome
#-------------------------------------------------------------------------------------
message("plotting coverage histogram")
map <- fread(
    cmd = paste("zcat", env$COVERAGE_FILE),
    col.names = c(
        "CHROM", # extract nodes fields
        "WINDOW",
        "COVERAGE"
    ),
    colClasses = c(
        "character", # PAF fields
        "integer",
        "integer"
    )
)
chroms <- map[, unique(CHROM)]
nChroms <- length(chroms)
chromIs <- as.list(1:nChroms)
names(chromIs) <- chroms
maxWindow <- map[, max(WINDOW)]
peakCoverage <- map[COVERAGE > 0, getmode(COVERAGE)]
ploidy <- 2
readsPerAllele <- peakCoverage / ploidy
maxY <- readsPerAllele * 5
# png(
#     # filename = paste(env$COVERAGE_FILE, "png", sep = "."),
#     filename = paste("/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svLong/HISTOGRAM", "png", sep = "."),
#     width = 3, 
#     height = 4, 
#     units = "in", 
#     pointsize = 8,
#     res = 600, 
#     type = "cairo"
# )
# dt <- map[, .(COUNT = .N), by = COVERAGE][order(COVERAGE)]
# plot(dt$COVERAGE, dt$COUNT, typ = "l", xlim = c(0, maxY))
# abline(v = readsPerAllele * 0:maxY)
# dev.off()
#=====================================================================================

#=====================================================================================
message("loading nodes")
nodes <- fread(
    env$NODES_FILE,
    col.names = c(
        "QNAME", # extract nodes fields
        "NODE1",
        "NODE2",
        "EDGE_TYPE",
        "MAPQ",
        "EVENT_SIZE",
        "INSERT_SIZE",
        "N_STRANDS"
    ),
    colClasses = c(
        "character", # PAF fields
        "integer",
        "integer",
        "character",
        "integer",
        "integer",
        "integer",
        "integer"
    )
)
str(nodes)

message("counting junctions")
junctions <- nodes[EDGE_TYPE != "A", .(
    EDGE_TYPE = EDGE_TYPE[1],
    N = .N,
    MAPQ = max(MAPQ),
    EVENT_SIZE = EVENT_SIZE[1]
), by = .(NODE1, NODE2)]
goodJunctions <- junctions[MAPQ > 55 & N > 3 & (EVENT_SIZE > 1e4 | EDGE_TYPE == "T")]
expanded1 <- parseSignedWindow(goodJunctions$NODE1)
setnames(expanded1, c("chromIndex1", "windowIndex1", "strand1"))
expanded2 <- parseSignedWindow(goodJunctions$NODE2)
setnames(expanded2, c("chromIndex2", "windowIndex2", "strand2"))
goodJunctions <- cbind(goodJunctions, expanded1, expanded2)
goodJunctions[, color := ifelse(EDGE_TYPE == "T", "red3", "orange2")]

str(junctions)
str(goodJunctions)
print(goodJunctions[, .(N = .N), by = "EDGE_TYPE"])


message("plotting genome")
png(
    # filename = paste(env$COVERAGE_FILE, "png", sep = "."),
    filename = paste("/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/glover/svPore/GENOME", "png", sep = "."),
    width = 12, 
    height = 0.5 * nChroms, 
    units = "in", 
    pointsize = 8,
    res = 300, 
    type = "cairo"
)

layout(matrix(1:nChroms, ncol = 1))
map[, {
    dt <- .SD[, .(
        x = mean(WINDOW),
        y = mean(COVERAGE)
    ), by = floor(WINDOW / 1e4)]
    par(mar = c(0.1, 4.1, 0.1, 0.1))
    plot(dt$x, dt$y, typ="p", pch = 19, cex = 0.25, col = rgb(0, 0, 1, 0.25),
         xaxt = "n",
         xlab = paste(CHROM, "Coordinate (bp)"), ylab = "Coverage", 
         xlim = c(1, maxWindow), ylim = c(0, maxY),
         yaxs = "i")
    abline(h = readsPerAllele * 0:maxY)
    j <- goodJunctions[chromIndex1 == chromIs[[CHROM]]]
    abline(v = j$windowIndex1 * 100, col = j$color)
    j <- goodJunctions[chromIndex2 == chromIs[[CHROM]]]
    abline(v = j$windowIndex2 * 100, col = j$color)
}, by = CHROM]
dev.off()
#=====================================================================================

#=====================================================================================
# message("collapsing nodes into molecules")
# molecules <- nodes[, .(
#     EDGES = paste(EDGE_TYPE, collapse = ""),       
#     N_EDGES = .N,
#     MIN_MAPQ = min(MAPQ),
#     MAX_MAPQ = max(MAPQ),
#     MIN_ALN_SIZE = min(EVENT_SIZE[EDGE_TYPE == "A"]),
#     MAX_ALN_SIZE = max(EVENT_SIZE[EDGE_TYPE == "A"]),    
#     MIN_SV_SIZE = if(.N == 1) 0 else min(EVENT_SIZE[EDGE_TYPE != "A"]),
#     MAX_SV_SIZE = if(.N == 1) 0 else max(EVENT_SIZE[EDGE_TYPE != "A"]),
#     N_STRANDS = N_STRANDS[1],
#     NODES = list(c(NODE1, NODE2[.N])),
#     PATH = getPathSignature(c(NODE1, NODE2[.N]))  ## IS THIS CORRECT?
 
# ), by = "QNAME"]

# str(nodes)
# str(molecules[N_EDGES > 1])
# str(molecules)

# # find duplex paths sequenced as two molecules (presumably rare)
# message("collapsing duplicates")
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

# str(molecules)

# # characterize individual SV junctions, count molecules

# # characterize composite SV paths, counts molecules

# # above required strand-aware comparison, similar to NODE_PATH1 == reverse(-NODE_PATH2)



#=====================================================================================
