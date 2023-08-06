#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)    
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
    flankSeq2 <- primers[2, substr(if(jxn$strand1 == "+") sequence else rcSequence, seqLength - flankSize2 + 1, seqLength)]
    paste0(flankSeq1, if(jxn$insertSize > 0) jxn$jxnSeq else "", flankSeq2)
}
junctions[, fakeSeq := unlist(mclapply(jxnKey, getJxnSeq, mc.cores = env$N_CPU))]
#=====================================================================================

# #########################
# dir <- "/nfs/turbo/umms-daviferg/jobFiles/TEW-svPore-test"
# eFile <- file.path(dir, "DEBUG.edges.rds")
# pFile <- file.path(dir, "DEBUG.primers.rds")
# jFile <- file.path(dir, "DEBUG.junctions.rds")
# nFile <- file.path(dir, "DEBUG.networks.rds")
# saveRDS(edges, file = eFile)
# saveRDS(primers, file = pFile)
# saveRDS(junctions, file = jFile)
# edges <- readRDS(eFile)
# primers <- readRDS(pFile)
# junctions <- readRDS(jFile)

#=====================================================================================
# extend junction matching by using edit distances between all unique junction pairs to create junction networks
#-------------------------------------------------------------------------------------
message("building junction networks from edit distances between unique junction sequences")
junctions <- junctions[order(-nMatchingSegments)]
junctions[, ":="(
    junctionI = 1:.N,
    parentJunctionI = NA_integer_
)]
isEditDistance1 <- function(qryJxnI, seedJxn){
    qryJxn <- junctions[qryJxnI]
    nodeDist <- sqrt((seedJxn$node1 - qryJxn$node1)**2 + (seedJxn$node2 - qryJxn$node2)**2) # first check if nodes are close enough to warrent edit distance calc
    if(nodeDist > 5) return(FALSE) # TODO: optimal threshold? expose as parameter?
    tryCatch({
        rawSeed <- charToRaw(seedJxn$fakeSeq) # adopt a simplified edit distance that builds on the specific construction of our sequences
        rawQry  <- charToRaw(qryJxn$fakeSeq)
        lengthSeed <- length(rawSeed)
        lengthQry  <- length(rawQry)
        minLength  <- min(lengthSeed, lengthQry)
        firstDiverged <- which(     rawSeed[1:minLength] !=      rawQry[1:minLength])[1]
        flank2        <- which(rev(rawSeed)[1:minLength] != rev(rawQry)[1:minLength])[1] - 1            
        adist(
            substr(seedJxn$fakeSeq, firstDiverged - 5, lengthSeed - flank2 + 5), # deliberately include extra shared bases if it helps adist find the best match
            substr(qryJxn$fakeSeq,  firstDiverged - 5, lengthQry  - flank2 + 5)
        ) == 1
    }, error = function(e){
        print(e)        
        message(paste(lengthSeed, lengthQry, minLength, firstDiverged, flank2))
        adist(seedJxn$fakeSeq, qryJxn$fakeSeq) == 1
    })
}
processSeedJxn <- function(seedJxnI, parentJxnI){
    junctions[seedJxnI, parentJunctionI := parentJxnI]
    qryJxnI <- junctions[is.na(parentJunctionI), junctionI]
    if(length(qryJxnI) == 0) return(NULL)

    # isEditDist1 <- unlist(lapply(qryJxnI, isEditDistance1, junctions[seedJxnI]))
    isEditDist1 <- unlist(mclapply(qryJxnI, isEditDistance1, junctions[seedJxnI], mc.cores = env$N_CPU))

    if(any(isEditDist1)){    
        junctions[qryJxnI[isEditDist1], parentJunctionI := parentJxnI]
        for(seedJxnI in qryJxnI[isEditDist1]) processSeedJxn(seedJxnI, parentJxnI) # iteratively use each matched child junction as a new seed in the network
    }
}
while(junctions[, any(is.na(parentJunctionI))]){
    parentJxnI <- junctions[is.na(parentJunctionI), junctionI[1]] # get the next network parent as the next most abundant but unassigned junction

    # parentJxnI <- 12

    nMatchingSegments <- junctions[parentJxnI, nMatchingSegments]
    if(nMatchingSegments == 1) break # stop assembling networks when down to singletons, they can't reasonably act as parents
    nRemaining <- junctions[, sum(is.na(parentJunctionI))]
    message(paste(parentJxnI, nMatchingSegments, nRemaining))
    processSeedJxn(parentJxnI, parentJxnI)
    ################
    # break
    # if(parentJxnI == 2) break
}
junctions[is.na(parentJunctionI), parentJunctionI := junctionI]
networks <- junctions[, .(
    networkKey = jxnKey[1], # first is the parent junction due to prior sorting
    nMatchingJunctions = .N,
    parentNMatchingSegments = nMatchingSegments[1],
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
    fakeSeq = fakeSeq[1]
), by = .(parentJunctionI)]
#=====================================================================================
 
#=====================================================================================
# report some useful tallies
#-------------------------------------------------------------------------------------
message("number of kept segments per read")
print(edges[, .(nKeptSegments = nKeptSegments[1]), by = .(qName)][, .N, by = .(nKeptSegments)][order(nKeptSegments)])
message("number of kept SV junctions per kept segment")
print(edges[, .(nKeptJunctions = nKeptJunctions[1]), by = .(segmentName)][, .N, by = .(nKeptJunctions)][order(nKeptJunctions)])
message("number of matching junctions per unique kept junction")
print(junctions[, .(nMatchingSegments), by = .(jxnKey)][, .N, by = .(nMatchingSegments)][order(nMatchingSegments)])
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
    networks, 
    paste(env$PARSE_PREFIX, "networks", "rds", sep = ".")
)
#=====================================================================================
