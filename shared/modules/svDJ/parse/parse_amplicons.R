#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)    
    library(bit64)  # provides support for 64-bit integers
    library(igraph) # provides support for directed graph analysis and clustering
    library(Rtsne)
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
        'PARSE_MODULE_DIR',
        'GENOME',
        'GENOME_FASTA',
        'SHM_DIR_WRK',
        'CHROM_FASTA_DIR',
        'DATA_NAME',
        'EXTRACT_PREFIX',
        'PARSE_PREFIX'
    ),
    integer = c(
        'N_CPU',
        'MIN_MAPQ',
        'MIN_ALIGNMENT_SIZE',
        'MIN_SV_SIZE',
        'JUNCTION_BANDWIDTH',
        'MIN_PARENT_COVERAGE',
        'MAX_EDIT_DISTANCE',
        'TSNE_PERPLEXITY'
    ), 
    double = c(
        'MIN_ALIGNMENT_IDENTITY',
        'MERGE_COVERAGE_SENSITIVITY',
        'TSNE_THETA'
    ),
    logical = c(
        'ADD_RARE_JUNCTIONS'
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general' , 'IUPAC', 'smith_waterman'))
sourceScripts(file.path(rUtilDir, 'genome'),   c('chroms','faidx'))
svPoreSharedDir <- file.path(env$MODULES_DIR, 'svPore')
sourceScripts(svPoreSharedDir, c(
    'constants','utilities','matching','filters'
))
sourceScripts(file.path(svPoreSharedDir, 'find'), c(
    'segments' 
))
sourceScripts(env$PARSE_MODULE_DIR, c(
    'filters'
))
setCanonicalChroms()
chromSizes <- loadChromSizes(env$WINDOW_SIZE)
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# load edges as prepared by `extract`
#-------------------------------------------------------------------------------------
message("loading edges")
edges <- readRDS(paste(env$ANALYZE_PREFIX, "edges", "rds", sep = "."))
#=====================================================================================

#=====================================================================================
# discover and order the amplicon primer node positions
#-------------------------------------------------------------------------------------
message("performing ab initio discovery of primer nodes")
reads <- edges[, .(outerNode1 = node1[1], outerNode2 = node2[.N]), by = .(readI)] # get the outer nodes, expecting two predominant strand orientations
primerNodes <- reads[, .(primerNode = abs(c(outerNode1, outerNode2)))][, .N, by = .(primerNode)][order(-N)] # set primer nodes as the most frequent, but distinct, outer nodes
primerNode1 <- primerNodes[1, primerNode]
i2 <- 2
primerNode2 <- primerNodes[i2, primerNode]
while(abs(primerNode2 - primerNode1) < 100){ # just in case the second most frequent primer node is a variant of primer1 (not likely...)
    i2 <- i2 + 1
    primerNode2 <- primerNodes[i2, primerNode]
}
primerNode1 <- reads[abs(outerNode1) == primerNode1, outerNode1[1]] # set the orientation of the primer nodes as the sign of the 5', i.e., first primer node
primerNode2 <- reads[abs(outerNode1) == primerNode2, outerNode1[1]]
swap <- if(sign(primerNode1) != sign(primerNode2)){
    sign(primerNode1) == -1 # order the primer nodes to hopefully (for ease of use) to yield a +/-, F/R pair
} else {
    abs(primerNode1) > abs(primerNode2) # otherwise order along genome so that primer order is deterministic
}
if(swap){
    tmp <- primerNode1
    primerNode1 <- primerNode2
    primerNode2 <- tmp
}
primers <- cbind(
    data.table(
        primer = 1:2,
        node = c(primerNode1, primerNode2)
    ),
    parseSignedNodes(chromSizes, c(primerNode1, primerNode2))
)
primers[, refPos := abs(refPos)]
print(primers)
rm(reads, primerNodes, swap)
#=====================================================================================

#=====================================================================================
# split reads at all chimeric junctions
# thus, each read might yield one or more non-chimeric segments
#-------------------------------------------------------------------------------------
message("matching junction nodes to primers")
isJunction <- getJunctionEdges(edges)
isNodeMatch <- function(node1, node2) abs(abs(node2) - abs(node1)) <= env$JUNCTION_BANDWIDTH 
edges[isJunction, ":="(
    node1IsPrimer = isNodeMatch(node1, primerNode1) | isNodeMatch(node1, primerNode2), # chimeric junctions could fuse primer1 to primer1, primer1 to primer2, etc.
    node2IsPrimer = isNodeMatch(node2, primerNode1) | isNodeMatch(node2, primerNode2)
)]
edges[isJunction, ":="(isChimeric = node1IsPrimer & node2IsPrimer)]
edges <- splitReadsToSegments(edges, getNonChimericJunctions, "splitting reads at chimeric junctions")
message(paste(edges[, length(unique(qName))], "reads ->", edges[, length(unique(segmentName))], "non-chimeric segments"))
rm(isJunction)
#=====================================================================================

#=====================================================================================
# purge segments that don't match primer1/primer2 at their outermost alignment nodes
#-------------------------------------------------------------------------------------
message("discarding segments that don't match primer1/primer2 amplicons")
segments <- edges[, .(outerNode1 = node1[1], outerNode2 = node2[.N]), by = .(segmentName)]
goodSegmentI <- segments[,
    (isNodeMatch(outerNode1, primerNode1) & isNodeMatch(outerNode2, primerNode2)) | 
    (isNodeMatch(outerNode1, primerNode2) & isNodeMatch(outerNode2, primerNode1))
]
edges <- edges[segmentName %in% segments[goodSegmentI, segmentName]]
message(paste(length(goodSegmentI), "segments ->", sum(goodSegmentI), "ampliconic, non-truncated segments"))
rm(segments, goodSegmentI)
#=====================================================================================

#=====================================================================================
# purge duplex segments from split chimeric molecules
# unlike svPore, do not purge duplex instances (yet) that arose from different reads
# as they may be truly independent sources molecules arising from biological or PCR replication
#-------------------------------------------------------------------------------------
message("discarding duplex segments (inverted segment orientations from the same read)")
edges[, nSegments := length(unique(segmentName)), by = .(qName)]
segmentPaths <- edges[nSegments > 1, .(
    nodePath    = paste(     c(node1, node2[.N]),  collapse = ":"),
    revNodePath = paste(-rev(c(node1, node2[.N])), collapse = ":")
), by = .(qName, segmentName)]
duplexSegmentsI <- segmentPaths[, c(
    sapply(1:(.N - 1), function(i1){
        any(sapply((i1 + 1):.N, function(i2) nodePath[i1] == revNodePath[i2])) # i.e., inverted path matches, here being strict about node positions
    }),
    FALSE # always keep the last segment, at least
), by = .(qName)][[2]]
duplexSegmentNames <- segmentPaths[duplexSegmentsI, segmentName]
edges <- edges[!(segmentName %in% duplexSegmentNames)]
message(paste(length(duplexSegmentsI), "segments from chimeric reads ->", sum(!duplexSegmentsI), "non-duplex segments"))
rm(segmentPaths, duplexSegmentsI, duplexSegmentNames)
#=====================================================================================

#=====================================================================================
# drop segments with no high-quality, matchable junctions
#-------------------------------------------------------------------------------------
message("dropping segments with no high-quality, usable junctions")
junctionI <- getJunctionEdges(edges)
keptJunctionI <- getMatchableJunctions(edges)
edges[, ":="( # set the count of total and kept junctions arising from each individual read segment
    nTotalJunctions = sum(junctionI[.I]),
    nKeptJunctions  = sum(keptJunctionI[.I])
), by = .(segmentName)]
segments <- edges[, .(kept = nKeptJunctions[1] > 0), by = .(segmentName)]
keptSegmentNames <- segments[kept == TRUE, unique(segmentName)]
edges <- edges[segmentName %in% keptSegmentNames]
edges[, nKeptSegments := length(unique(segmentName)), by = .(qName)]
message(paste(nrow(segments), "segments ->", segments[, sum(kept)], "kept segments"))
rm(junctionI, keptJunctionI, segments, keptSegmentNames)
#=====================================================================================

#=====================================================================================
# perform _exact_ junction matching, i.e., any base difference at a junction will lead to match failure
#-------------------------------------------------------------------------------------
message("strictly grouping junctions based on exact junction sequence matches")
edges[, isCanonical := isNodeMatch(node1[1], primerNode1), by = .(segmentName)] # canonical segments are primer1/primer2, non-canonical are primer2/primer1
isSingletonMatchable <- getSingletonMatchableJunctions(edges)
edges[isSingletonMatchable, jxnKey := ifelse(
    isCanonical, # microhomology bases are NOT used in jxnKey - the aligner already decided they were a match to reference
    paste( node1,  node2, insertSize, ifelse(insertSize > 0,    jxnSeq,  ""), sep = ":"),
    paste(-node2, -node1, insertSize, ifelse(insertSize > 0, rc(jxnSeq), ""), sep = ":")
)]
edges[isSingletonMatchable, cJxnSeq := ifelse(isCanonical, jxnSeq, rc(jxnSeq))] # keep information that may not be present in jxnKey
junctions <- edges[isSingletonMatchable, { # restrict to single-junction segments to facilitate edit distance, below
    x <- strsplit(jxnKey, ":")[[1]]
    .(
        nMatchingSegments = .N, # valid, since each segment has a single junction
        nCanonical = sum(isCanonical),
        nNonCanonical = sum(!isCanonical),
        node1 = as.integer64(x[1]),
        node2 = as.integer64(x[2]),
        insertSize = as.integer(x[3]),
        jxnSeq = paste(unique(cJxnSeq), collapse = ",") # preserves any variability found in microhomologous bases omitted from jxnKey
    )                                                   # always a single value for insertions since those jxnBases are part of jxnKey
}, by = .(jxnKey)]
junctions <- cbind(
    junctions, 
    junctions[, parseSignedNodes(chromSizes, node1, 1)], 
    junctions[, parseSignedNodes(chromSizes, node2, 2)]
)
rm(isSingletonMatchable)
#=====================================================================================

#=====================================================================================
# prepare for junction network analysis
#-------------------------------------------------------------------------------------
message("extracting genomic sequences at primer nodes")
primer1RefPos <- primers[1, refPos]
primer2RefPos <- primers[2, refPos]
maxFlank1 <- junctions[, max(abs(refPos1 - primer1RefPos) + 1)]
maxFlank2 <- junctions[, max(abs(refPos2 - primer2RefPos) + 1)]
loadFaidx()
primers[, sequence := {
    maxSize <- if(primer == 1) maxFlank1 else maxFlank2
    range <- if(strand == "+") c(refPos, refPos + maxSize - 1) 
                          else c(refPos - maxSize + 1, refPos)
    getRefSeq(chrom, range)
}, by = .(primer)]
primers[, rcSequence := rc(sequence)]
primers[, seqLength := nchar(sequence)]
unlink(genome_index$chromFile)
rm(genome_index, primer1RefPos, primer2RefPos, maxFlank1, maxFlank2)

# construct a fictitious, idealized read sequence in which any bases aligned as reference by the aligner
# are returned as reference (not read) bases, to ignore inconsequential errors in flanking sequences or microhomologies
# in contrast, inserted, non-reference bases are included as found in the read; we must still asses the impact of errors at these junction bases
message("assembling idealized junction sequences for junction network analysis")
getJxnSeq <- function(jxnKey_){
    jxn <- junctions[jxnKey == jxnKey_]
    microhomologyAdjustment <- max(0, -jxn$insertSize)
    flankSize1 <- abs(jxn$refPos1 - primers[1, refPos]) + 1 - microhomologyAdjustment
    flankSize2 <- abs(jxn$refPos2 - primers[2, refPos]) + 1 
    flankSeq1 <- primers[1, substr(if(jxn$strand1 == "+") sequence else rcSequence, 1, flankSize1)]
    flankSeq2 <- primers[2, substr(if(jxn$strand2 == "+") sequence else rcSequence, seqLength - flankSize2 + 1, seqLength)]
    paste0(flankSeq1, if(jxn$insertSize > 0) jxn$jxnSeq else "", flankSeq2)
}
junctions[, fakeSeq := unlist(mclapply(jxnKey, getJxnSeq, mc.cores = env$N_CPU))]
junctions[, fakeLength := nchar(fakeSeq)]

message("sorting junctions by number of matching read segments")
junctions <- junctions[order(-nMatchingSegments, -fakeLength)]
junctions[, junctionI := 1:.N]
#=====================================================================================

#=====================================================================================
# construct and analyze a directed graph where all junctions above a coverage threshold are nodes
# and edges are the edit distance between them, pointing from lower to higher coverage junctions
#-------------------------------------------------------------------------------------
message("calculating pairwise edit distances between unique junction sequences")
parentCandidates <- junctions[nMatchingSegments >= env$MIN_PARENT_COVERAGE]
nParentCandidates <- nrow(parentCandidates)
nParentNonCandidates <- nrow(junctions) - nParentCandidates
adistFlankBases <- 10
junctionEditDistance <- function(qry, ref, ...){ 
    if(qry == ref) return(0) # unclear why this happens, presumably minimap2 isn't determinsitic in its clipping
    qryN <- nchar(qry)
    refN <- nchar(ref)
    deltaN <- abs(qryN - refN)
    if(deltaN > env$MAX_EDIT_DISTANCE) return(deltaN) # shortcut when gaps are expected to predominate the edit distance
    qryRaw <- charToRaw(qry) # speed up adist by trimming shared flanks
    refRaw <- charToRaw(ref)
    minN <- min(qryN, refN)
    firstDiverged <- which(     qryRaw[1:minN] !=      refRaw[1:minN])[1]
    flank2        <- which(rev(qryRaw)[1:minN] != rev(refRaw)[1:minN])[1] - 1 
    qry <- substr(qry, firstDiverged - adistFlankBases, qryN - flank2 + adistFlankBases) # keep small anchor flanks to ensure edits are in the junction
    ref <- substr(ref, firstDiverged - adistFlankBases, refN - flank2 + adistFlankBases)
    adist(qry, ref)
}
pDist <- function(x, distFn, diag = FALSE, upper = FALSE){ # parallelized dist function   TODO: move this to a utilities script
    n <- length(x)
    v <- unlist(mclapply(1:(n-1), function(j) {
        i <- (j + 1):n
        sapply(x[i], distFn, x[j], i, j)
    }, mc.cores = env$N_CPU))
    attr(v, "Size") <- n 
    attr(v, "Diag") <- diag
    attr(v, "Upper") <- upper
    class(v) <- "dist"
    names(v) <- NULL
    v
}
candidateDistances <- unname(pDist(parentCandidates$fakeSeq, junctionEditDistance))

message("calculating T-SNE on the distance matrix")
tsne <- Rtsne(
    candidateDistances,
    perplexity = env$TSNE_PERPLEXITY,
    theta = env$TSNE_THETA,
    check_duplicates = FALSE,
    pca = FALSE,
    is_distance = TRUE,
    num_threads = env$N_CPU
)

message("building and analyzing the hierarchical junction network")
distanceMatrix <- as.matrix(candidateDistances)
parentJunctionI  <- integer()
getNetworkTable <- function(i, j, distance){
    if(i == j) parentJunctionI <<- c(parentJunctionI, j)
    data.table(junctionI = i, parentJunctionI = j, parentDistance = distance)
}
setBestParent <- function(i, parentDistances){
    bestParentI <- which(parentDistances == min(parentDistances))[1]
    getNetworkTable(i, parentJunctionI[bestParentI], parentDistances[bestParentI])     
}

message("  identifying parent junctions of junction communities")
membership <- do.call(rbind, lapply(1:nParentCandidates, function(i){
    if(i == 1) getNetworkTable(i, i, 0) else { # the most frequent junction sequence is always an parent
        parentDistances <- distanceMatrix[i, parentJunctionI]
        parentSizes <- parentCandidates[parentJunctionI, nMatchingSegments]    
        querySize <- parentCandidates[i, nMatchingSegments]
        relativeSizeAdjustment <- (parentSizes - querySize) / parentSizes # require greater distance for similar sized parent candidates
        relativeSizeAdjustment <- 1 - (1 - relativeSizeAdjustment) * env$MERGE_COVERAGE_SENSITIVITY # allow user to adjust the merge sensitivity
        parentMatches <- parentDistances <= env$MAX_EDIT_DISTANCE * relativeSizeAdjustment
        if(sum(parentMatches) == 0) getNetworkTable(i, i, 0) # create a new parent
        else setBestParent(i, parentDistances)   
    }
}))
if(nParentNonCandidates > 0) {
    if(env$ADD_RARE_JUNCTIONS) message("  adding evidence from low-coverage junctions")
    nonCandidateI <- (nParentCandidates + 1):nrow(junctions)
    membership <- rbind(
        membership, 
        if(env$ADD_RARE_JUNCTIONS) do.call(rbind, mclapply(nonCandidateI, function(i){
            parentDistances <- sapply(parentJunctionI, function(j) junctionEditDistance(junctions[i, fakeSeq], junctions[j, fakeSeq]))
            parentMatches <- parentDistances <= env$MAX_EDIT_DISTANCE
            if(sum(parentMatches) == 0) data.table(junctionI = i, parentJunctionI = NA, parentDistance = NA) # drop/mask low-coverage junctions with no parent match
            else setBestParent(i, parentDistances)   
        }, mc.cores = env$N_CPU)) else data.table(junctionI = nonCandidateI, parentJunctionI = NA, parentDistance = NA)
    )
}

message("  making network summary plot")
nParentJunctions <- length(parentJunctionI)
cols <- sample(rainbow(nParentJunctions))
png(filename = paste(env$PARSE_PREFIX, "junction_graph.png",  sep = "."), 
    width = 3, height = 3, units = "in", pointsize = 8,
    res = 600, type = "cairo")
plot(tsne$Y, pch = 1, cex = parentCandidates$nMatchingSegments / max(parentCandidates$nMatchingSegments) * 3, 
     col = ifelse(membership[, junctionI == parentJunctionI], "black", cols[membership$parentJunctionI]),
     main = env$DATA_NAME)
invisible(dev.off())

message("  aligning each junction to its network parent")
junctions <- merge(junctions, membership, by = "junctionI", all.x = TRUE)
maxShift <- env$MAX_EDIT_DISTANCE * 2
for(parentJxnI in parentJunctionI){
    networkJxnI <- junctions[parentJunctionI == parentJxnI, junctionI]
    junctions[networkJxnI, junctionOnParent := unlist(mclapply(networkJxnI, function(jxnI){
        paste0(smith_waterman(junctions[jxnI, fakeSeq], junctions[parentJxnI, fakeSeq], fast = TRUE)$qryOnRef, collapse = "")
    }, mc.cores = env$N_CPU))]
}

message("  assembling the networks table")
networks <- junctions[!is.na(parentJunctionI)][, .(
    networkKey = jxnKey[1], # first junction is the parent junction due to prior sorting
    maxDistance = max(parentDistance, na.rm = TRUE),
    nMatchingJunctions = .N,
    parentNMatchingSegments = nMatchingSegments[1],
    nextNMatchingSegments = if(.N == 1) NA_integer_ else nMatchingSegments[2],
    nMatchingSegments = sum(nMatchingSegments),
    nCanonical = sum(nCanonical),
    nNonCanonical = sum(nNonCanonical),
    node1 = node1[1],
    node2 = node2[1],
    insertSize = insertSize[1],
    jxnSeq = jxnSeq[1],
    chrom1 = chrom1[1],
    chromIndex1 = chromIndex1[1],
    refPos1 = refPos1[1],
    strand1 = strand1[1],
    chrom2 = chrom2[1],    
    chromIndex2 = chromIndex2[1],
    refPos2 = refPos2[1],
    strand2 = strand2[1],
    fakeSeq = fakeSeq[1],
    fakeLength = fakeLength[1]
), by = .(parentJunctionI)]
#=====================================================================================

#=====================================================================================
# report some useful tallies
#-------------------------------------------------------------------------------------
message("number of kept segments per read")
print(edges[, .(nKeptSegments = nKeptSegments[1]), by = .(qName)][, .N, by = .(nKeptSegments)][order(nKeptSegments)])
message("number of kept SV junctions per kept segment")
print(edges[, .(nKeptJunctions = nKeptJunctions[1]), by = .(segmentName)][, .N, by = .(nKeptJunctions)][order(nKeptJunctions)])
message("number of matching segments per unique kept junction")
print(junctions[, .N, by = .(nMatchingSegments)][order(nMatchingSegments)])
message("number of matching segments per junction network")
print(networks[, .N, by = .(nMatchingSegments)][order(nMatchingSegments)])
message("number of junctions per junction network")
print(networks[, .N, by = .(nMatchingJunctions)][order(nMatchingJunctions)])
message("maximum edit distance per junction network")
print(networks[, .N, by = .(maxDistance)][order(maxDistance)])
#=====================================================================================

#=====================================================================================
# make a simple sample manifest and commit results to disk
#-------------------------------------------------------------------------------------
message("creating sample manifest")    
nEdges <- nrow(edges)
nJunctions <- nrow(junctions)
nNetworks <- nrow(networks)
message(paste("nEdges", nEdges))
message(paste("nJunctions", nJunctions))
message(paste("nNetworks", nNetworks))
isSingletonMatchable <- getSingletonMatchableJunctions(edges)
manifest <- edges[isSingletonMatchable, .(
    Project = env$DATA_NAME,
    Sample_ID = sample[1],
    Description = sample[1],
    Yield = .N,
    nSingletonMatchable = .N,
    nEdges = nEdges,
    nJunctions = nJunctions,
    nNetworks = nNetworks
)]
fwrite(
    manifest,
    file = paste(env$PARSE_PREFIX, 'sample_manifest', 'csv', sep = "."), 
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
    paste(env$PARSE_PREFIX, "chromosome_metadata", "rds", sep = ".")
)

# save results for the app
message("saving parse results") 
saveRDS(
    primers, 
    paste(env$PARSE_PREFIX, "primers", "rds", sep = ".")
)
saveRDS(
    edges, 
    paste(env$PARSE_PREFIX, "edges", "rds", sep = ".")
)
saveRDS(
    junctions, 
    paste(env$PARSE_PREFIX, "junctions", "rds", sep = ".")
)
saveRDS(
    distances, 
    paste(env$PARSE_PREFIX, "distances", "rds", sep = ".")
)
saveRDS(
    networks, 
    paste(env$PARSE_PREFIX, "networks", "rds", sep = ".")
)
#=====================================================================================


# # EARLY CODE USING UMI-DERIVED DIRECTIONAL NETWORK APPROACH
# # could probably still do it this way if size thresholding was introduced 
# # but approach above is more efficient and accurate
# junctions[, ":="(
#     junctionI = 1:.N,
#     parentJunctionI = NA_integer_,
#     parentEditDistance = 0
# )]
# isEditDistance1 <- function(qryJxnI, seedJxn){
#     qryJxn <- junctions[qryJxnI]
#     nodeDist <- sqrt((seedJxn$refPos1 - qryJxn$refPos1)**2 + (seedJxn$refPos2 - qryJxn$refPos2)**2) # first check if nodes are close enough to warrent edit distance calc
#     if(nodeDist > env$JUNCTION_BANDWIDTH || abs(qryJxn$fakeLength - seedJxn$fakeLength) > 2) return(FALSE)
#     rawSeed <- charToRaw(seedJxn$fakeSeq) # adopt a simplified edit distance that builds on the specific construction of our sequences
#     rawQry  <- charToRaw(qryJxn$fakeSeq)
#     lengthSeed <- length(rawSeed)
#     lengthQry  <- length(rawQry)
#     minLength  <- min(lengthSeed, lengthQry)
#     firstDiverged <- which(     rawSeed[1:minLength] !=      rawQry[1:minLength])[1]
#     flank2        <- which(rev(rawSeed)[1:minLength] != rev(rawQry)[1:minLength])[1] - 1 
#     editDistance <- adist(
#         substr(seedJxn$fakeSeq, firstDiverged - 5, lengthSeed - flank2 + 5), # deliberately include extra shared bases to help adist find the best edit path
#         substr(qryJxn$fakeSeq,  firstDiverged - 5, lengthQry  - flank2 + 5)
#     )
#     if(is.na(editDistance)) editDistance <- adist(seedJxn$fakeSeq, qryJxn$fakeSeq) # catch rare problems with vector overruns in assembly above
#     editDistance <= 1 # should never be zero...
# }
# processSeedJxn <- function(seedJxnI, parentJxnI, parentEditDistance_){
#     seedJxn <- junctions[seedJxnI]    
#     qryJxnIs <- junctions[is.na(parentJunctionI), junctionI]
#     if(length(qryJxnIs) == 0) return(NULL)
#     isEditDist1 <- unlist(mclapply(qryJxnIs, isEditDistance1, seedJxn, mc.cores = env$N_CPU))
#     if(any(isEditDist1)){  
#         junctions[qryJxnIs[isEditDist1], ":="(parentJunctionI = parentJxnI, parentEditDistance = parentEditDistance_)]
#         for(childJxnI in qryJxnIs[isEditDist1]) processSeedJxn(childJxnI, parentJxnI, parentEditDistance_ + 1) # iteratively use each matched child junction as a new seed in the network
#     }
#     NULL
# }
# message(paste("parentJxnI", "nMatchingSegments", "nRemaining"))
# pendingJxns <- junctions[, is.na(parentJunctionI)]
# maxShift <- env$JUNCTION_BANDWIDTH
# while(any(pendingJxns)){
#     parentJxnI <- which(pendingJxns)[1] # get the next network parent as the next most abundant but unassigned junction   
#     nMatchingSegments <- junctions[parentJxnI, nMatchingSegments]
#     if(nMatchingSegments == 1) break # stop assembling networks when down to singletons, they can't reasonably act as parents
#     message(paste(parentJxnI, nMatchingSegments, sum(pendingJxns)))
#     junctions[parentJxnI, ":="(parentJunctionI = parentJxnI, parentEditDistance = 0)]
#     processSeedJxn(parentJxnI, parentJxnI, 1)
#     pendingJxns <- junctions[, is.na(parentJunctionI)]
# }
# junctions[pendingJxns, ":="(parentJunctionI = junctionI)]

# # EARLY CODE USING Smith-Waterman distance
# # abandone because eventually was simply turning SW into adist, just use adist
# matchScore          <-  1
# mismatchPenalty     <- -1 # -1.5
# gapOpenPenalty      <- -1.001 # -2.501 # 0.001 ajustment gives slight preference to not opening a single-base terminal gap
# gapExtensionPenalty <- -1
# maxShift            <- env$JUNCTION_BANDWIDTH
# pairedBaseScores <- initializePairScores()
# swFlankBases <- 10
# identityDistance <- 0.01 # not zero, since that causes igraph to not create the edge
# SWDist <- function(qry, ref, qryI, refI){ # 2/1 order reflect the call from pDist
#     if(qry == ref) return(identityDistance) # unclear why this happens, presumably minimap2 isn't determinsitic in its clipping, etc.
#     qryN <- nchar(qry)
#     refN <- nchar(ref)
#     deltaN <- abs(qryN - refN)
#     if(deltaN > maxShift) return(abs( # fast approximation when we expect gaps to predominate the score
#         gapOpenPenalty * (deltaN / maxShift) + # scale gap open penalties to a fraction of total gap bases
#         deltaN * gapExtensionPenalty
#     )) 
#     qryRaw <- charToRaw(qry) # speed up fast Smith-Waterman even more by trimming shared flanks
#     refRaw <- charToRaw(ref)
#     minN <- min(qryN, refN)
#     firstDiverged <- which(     qryRaw[1:minN] !=      refRaw[1:minN])[1]
#     flank2        <- which(rev(qryRaw)[1:minN] != rev(refRaw)[1:minN])[1] - 1 
#     qry <- substr(qry, firstDiverged - swFlankBases, qryN - flank2 + swFlankBases)
#     ref <- substr(ref, firstDiverged - swFlankBases, refN - flank2 + swFlankBases)
#     max(nchar(qry), nchar(ref)) - smith_waterman(qry, ref, fast = TRUE)$bestScore # thus, the net penalty vs. the best possible score
# }
# distances <- unname(pDist(junctions$fakeSeq, SWDist))

# # EARLY CODE USING IGRAPH
# # abandoned because no community-finding methods use directed networks, edge weights and node sizes
# message("constructing an igraph from the distance matrix")
# swGraphThreshold <- 6
# distanceMatrix <- as.matrix(distances)
# dm <- distanceMatrix
# dm[upper.tri(dm)] <- 0 # force the graph to be directional toward the most highly represented junction sequence nodes
# diag(dm) <- 0
# dm[dm > swGraphThreshold] <- 0
# junctionGraph <- graph_from_adjacency_matrix(dm, mode = "directed", weighted = TRUE, diag = FALSE)

# message("assigning junction communities")
# junctionI <- 1:nrow(junctions)
# ancestorJunctions <- lapply(junctionI, function(i) subcomponent(junctionGraph, i, mode = "out"))
# isParentJunction <- sapply(ancestorJunctions, length) == 1
# isParentJunction <- sapply(junctionI, function(i) isParentJunction[i] & junctions[i, nMatchingSegments > 5])
# nParentJunctions <- sum(isParentJunction) 
# membership <- sapply(junctionI, function(i){
#     if(isParentJunction[i]) i
#     else {
#         v <- ancestorJunctions[[i]]
#         parentJunctions <- v[isParentJunction[v]]
#         if(length(parentJunctions) == 0) i
#         else if(length(parentJunctions) == 1) parentJunctions
#         else {
#             parentDistances <- sapply(parentJunctions, function(j) as.vector(distanceMatrix[i, j]))
#             parentJunctions[which(parentDistances == min(parentDistances))[1]]
#         }
#     }
# })
