# sort and group molecules based on the SVs they carry

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("  initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOMEX_MODULES_DIR',
        'OUTPUT_DIR',
        'DATA_NAME',
        'DATA_FILE_PREFIX',
        'INTERIM_FILE',
        'MOLECULE_TYPES_FILE',
        'JUNCTIONS_FILE',
        'AMPLICONS_FILE'
    ),
    integer = c(
        'N_CPU',
        'MIN_SV_SIZE'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'), c('general','IUPAC','smith_waterman'))
#-------------------------------------------------------------------------------------
# set some constants
edgesColumns <- c(
    "molId", # was called molId upstream, will be become the type id once an index is selected
    "chrom1", 
    "strand1",
    "pos1",
    "seq1",
    "qual1",    
    "readN1",
    "chrom2",
    "strand2",
    "pos2",
    "seq2",
    "qual2",    
    "readN2",
    "mapQ",
    "cigar",
    "alnQ",
    "edgeType",
    "eventSize",
    "insertSize",
    "jxnBases",
    "ampliconId", 
    "mergeLevel", 
    "overlap", 
    "isReference", 
    "nReadPairs",
    "baseQual"
)
edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    MERGE_FAILURE = "M",
    PROPER        = "P"
)
inDelTypes <- c("D","I")
#=====================================================================================

#=====================================================================================
message("  loading edges")
edges <- fread(env$INTERIM_FILE, sep = "\t", header = FALSE)
setnames(edges, edgesColumns)
message("  calculating derivative values")
isJunction      <- edges[, edgeType != edgeTypes$ALIGNMENT]
isMergeFailure  <- edges[, edgeType == edgeTypes$MERGE_FAILURE]
isTooSmall <- edges[, 
    edgeType %in% inDelTypes & # true junctions are masked if inline and too small
    eventSize  < env$MIN_SV_SIZE & 
    insertSize < env$MIN_SV_SIZE 
]
isRejectedIndel <- edges[, # called bandwidth in other svX pipelines, the net reference offset of a paired indel
    edgeType == edgeTypes$DELETION & # per extract.sh, edges are called deletion/D unless eventSize == 0
    abs(eventSize - insertSize) < env$MIN_SV_SIZE
]
edges[, ":="(molKey = paste(ampliconId, molId, sep = ":"))] # used for simplified sorting and grouping
isFailedBandwidth <- edges[, # bandwidth failure fails all junctions in the molecule
    rep(if(mergeLevel[1] > 0){
        refSpan <- abs(pos1[1] - pos2[.N])
        qrySpan <- nchar(seq1[1])
        abs(refSpan - qrySpan) < env$MIN_SV_SIZE
    } else FALSE, .N), # cannot meaningfully check bandwith if reads did not merge
    by = .(molKey)
][[2]]
edges[, ":="(
    i = 1:.N,
    mergeFailure = isMergeFailure,
    tooSmall = isTooSmall,
    rejectedIndel = isRejectedIndel,
    failedBandwidth = isFailedBandwidth,
    maskedJunction = isMergeFailure | isTooSmall | isRejectedIndel | isFailedBandwidth,
    candidateJunction = isJunction & !isMergeFailure & !isTooSmall & !isRejectedIndel & !isFailedBandwidth,
    node1 = paste(chrom1, strand1, pos1, sep = "/"),
    node2 = paste(chrom2, strand2, pos2, sep = "/")
)]
bQ <- edges$baseQual
edges[isJunction, ":="( # junction qualities are the lower of the two flanking alignment qualities
    baseQual = pmin(bQ[i - 1], bQ[i + 1]), # mapQ and alnQ were already set properly upstream
    insertBases = ifelse(insertSize > 0, jxnBases, "*") # insert (but not microhomology) bases for junction grouping    
)]
edges[isMergeFailure, ":="( # estimate overlap of merge failures
    overlap = ifelse(chrom1 == chrom2, pos1 - pos2 + 1, NA) # thus, a positive overlap for the primary type of expected merge failure
)]
setkey(edges, molKey, i)
rm(isJunction, isMergeFailure, isRejectedIndel, bQ)
#=====================================================================================

#=====================================================================================
message("  counting unique molecules")
nMolecules     <- edges[,                       length(unique(molKey))]
nMergeFailure  <- edges[mergeFailure == TRUE,   length(unique(molKey))]
nMasked        <- edges[maskedJunction == TRUE, length(unique(molKey))]
message(paste0("    ", paste(nMolecules,     "unique molecule sequences",     sep = "\t")))
message(paste0("    ", paste(nMergeFailure,  "molecules had merge failures",  sep = "\t")))
message(paste0("    ", paste(nMasked,        "molecules a masked junction",  sep = "\t")))
#=====================================================================================

#=====================================================================================
message("  collapsing edges into molecules")
molecules <- edges[, {
    isAlignment <- edgeType == edgeTypes$ALIGNMENT    
    isMerged <- mergeLevel[1] > 0
    interleave_ <- function(value1, value2, filter){ # fields that may have different values on the two sides of a junction
        paste(c(rbind(value1[filter], value2[filter])), collapse = ":")
    }
    paste_ <- function(value, filter){ # fields that have the same value for the entire junction
        paste(value[filter], collapse = ":")
    }
    min_ <- function(value, filter){
        min(value[filter], na.rm = TRUE)
    }
    .(
        # molecule-level information, the same for every edge
        # molKey
        ampliconId = ampliconId[1], 
        molId = molId[1],  
        isReference =  as.logical(isReference[1]),    
        nReadPairs = nReadPairs[1],    
        mergeLevel = mergeLevel[1], 
        overlap = overlap[1],

        # path-level information, one concatenated/aggregated value per molecule over all unmasked junctions
        # for amplicon sequencing, the grouping path need not include the outer alignment nodes        
        path = interleave_(node1, node2,     candidateJunction),
        pathEdgeTypes   = paste_(edgeType,   candidateJunction),
        pathInsertSizes = paste_(insertSize, candidateJunction),
        pathJxnBases    = paste_(jxnBases,   candidateJunction),
        pathInsertBases = paste_(insertBases,candidateJunction),
        minMapQ     = min_(mapQ,     isAlignment), # how bad the worst part of a path was
        minAlnQ     = min_(alnQ,     isAlignment),
        minBaseQual = min_(baseQual, isAlignment),  

        # node-level information, one value per unique node, including the outermost nodes
        poss = paste(c(pos1, pos2[.N]), collapse = ":"),

        # junction-level information, concatenated across all junctions, one value per junction
        edgeTypes   = paste_(edgeType,   !isAlignment), 
        eventSizes  = paste_(eventSize,  !isAlignment),
        insertSizes = paste_(insertSize, !isAlignment),
        jxnBases    = paste_(jxnBases,   !isAlignment),
        insertBases = paste_(insertBases,!isAlignment),
        mapQs       = paste_(mapQ,       !isAlignment), # quality values already set as the minimum of the two flanks
        alnQs       = paste_(alnQ,       !isAlignment),
        baseQuals   = paste_(baseQual,   !isAlignment),

        # alignment-level information, one value per alignment, all values here match between node1 and node2 of alignments
        chroms      = paste_(chrom1,  isAlignment), 
        strands     = paste_(strand1, isAlignment),
        cigars      = paste_(cigar,   isAlignment),         
        readNs      = paste_(readN1,  isAlignment),

        # read-level information
        seq1  = seq1[1], 
        seq2  = if(isMerged) "*" else seq2[.N], # there are only ever one merged or two unmerged reads
        qual1 = qual1[1],
        qual2 = if(isMerged) "*" else qual2[.N]
    )
}, by = .(molKey)]
#=====================================================================================

#=====================================================================================
message("  aggregating molecules into types, i.e., shared SV paths")
moleculeTypes <- molecules[
    order(
        ampliconId, path, pathInsertSizes, pathInsertBases, # first columns are for moleculeType grouping
        -isReference, -mergeLevel, -nReadPairs, -minBaseQual, -minMapQ # last columns are   
    )
][, .(
        # molecule-level information
        # ampliconId 
        molKey      = molKey[1],  # one molecule id selected as the index for the type
        molId       = molId[1],        
        isReference = isReference[1],
        nReadPairs  = sum(nReadPairs), # continue summing all original _read pairs_ that has this type 
        nMolecules  = .N, # the number of distinct read pair _sequences_ that had this type  
        mergeLevel  = mergeLevel[1],
        overlap     = na.omit(overlap)[1],

        # path-level information
        # path
        # pathInsertSizes
        # pathInsertBases
        pathEdgeTypes = pathEdgeTypes[1],
        pathJxnBases  = pathJxnBases[1],
        minMapQ     = max(minMapQ), # thus, take the worst segment score per molecule, but the best of those values across all molecules matching a type
        minAlnQ     = max(minAlnQ),
        minBaseQual = max(minBaseQual),

        # node-level information
        poss        = poss[1],     

        # junction-level information
        edgeTypes   = edgeTypes[1],
        eventSizes  = eventSizes[1],   
        insertSizes = insertSizes[1],  
        jxnBases    = jxnBases[1],   
        insertBases = insertBases[1],  
        mapQs       = mapQs[1],
        alnQs       = alnQs[1],
        baseQuals   = baseQuals[1],

        # alignment-level information
        chroms      = chroms[1],
        strands     = strands[1],
        cigars      = cigars[1],
        readNs      = readNs[1],         

        # read-level information
        seq1        = seq1[1],
        seq2        = seq2[1],        
        qual1       = qual1[1],
        qual2       = qual2[1]
), by = .(ampliconId, path, pathInsertSizes, pathInsertBases)][
    order(ampliconId, -nReadPairs) # sort by molecule type frequency
]
#=====================================================================================

#=====================================================================================
message("  tabulating unique SV junctions found in one or more molecules")
molTypeKeys <- moleculeTypes[, molKey]
edges[, isIndexMol := molKey %in% molTypeKeys]
junctions <- edges[edgeType != edgeTypes$ALIGNMENT & maskedJunction == FALSE, .( # only tabulate candidate junction, not merge failures our under-sized SVs
    # ampliconId
    # node1
    # node2
    # insertSize
    # insertBases
    nReadPairs  = sum(nReadPairs),
    nMolecules  = length(unique(molKey)), # beware of possibility of rare molecules with a recurring junction
    nMolTypes   = length(unique(molKey[isIndexMol])), 
    chrom1      = chrom1[1],
    strand1     = strand1[1],
    pos1        = pos1[1],
    chrom2      = chrom2[1],
    strand2     = strand2[1],
    pos2        = pos2[1],
    mapQ        = max(mapQ),
    alnQ        = max(alnQ),
    baseQual    = max(baseQual),    
    edgeType    = edgeType[1],
    eventSize   = eventSize[1],
    jxnBases    = jxnBases[1],
    mergeLevel  = max(mergeLevel),
    molTypeKeys = list(unique(molKey[isIndexMol]))
), by = .(ampliconId, node1, node2, insertSize, insertBases)] # include insertSize and insertBases in junction definitons
#=====================================================================================

#=====================================================================================
message("  saving outputs")
saveRDS(moleculeTypes, file = env$MOLECULE_TYPES_FILE)
saveRDS(junctions, file = env$JUNCTIONS_FILE)
write.table(
    data.table(
        Project     = basename(env$OUTPUT_DIR),
        Sample_ID   = env$DATA_NAME,
        Description = env$DATA_NAME
    ),
    file = env$MANIFEST_FILE, 
    quote = TRUE, 
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)
#=====================================================================================

# #=====================================================================================
# # append the original raw unmerged reads for the representative of each molecule type
# message("  extracting representative original unmerged reads")
# tabixFile <- paste0(env$DATA_FILE_PREFIX, ".indexed_reads.bgz")
# tbx <- Rsamtools::TabixFile(tabixFile)
# reads <- do.call(rbind, mclapply(moleculeTypes[, molId], function(molId){
#     index <- floor(molId / 100) + 1
#     i <- molId - (index - 1) * 100 + 1
#     x <- Rsamtools::scanTabix(tbx, param = GenomicRanges::GRanges("X", IRanges::IRanges(index, width = 1)))
#     x <- strsplit(x[[1]][i], "\t")[[1]]
#     data.table(
#         read1Seq  = x[4],
#         read2Seq  = x[6]
#     )
# }, mc.cores = env$N_CPU))
# moleculeTypes <- cbind(moleculeTypes, reads)
# moleculeTypes[, ":="( # which read matches seq1 (which is always the left end)
#     strand = as.integer(
#         as.vector(adist(substr(read1Seq, 1, 20), substr(seq1, 1, 20))) >=
#         as.vector(adist(substr(read2Seq, 1, 20), substr(seq1, 1, 20)))
#     )
# ), by = "molKey"]
# message("  re-aligning original reads to the amplicon reference")
# amplicons <- fread(env$AMPLICONS_FILE, col.names = c(
#     "amplicon","type","nReadPairs",
#     "chrom1","side1","pos1","ref1","primer1",
#     "chrom2","side2","pos2","ref2","primer2"
# ))
# moleculeTypes <- merge(
#     moleculeTypes, 
#     do.call(rbind, mclapply(1:nrow(moleculeTypes), function(i){
#         moleculeTypes[i, {
#             amp <- amplicons[amplicon == ampliconId]
#             swL <- smith_waterman(
#                 if(strand == 0) read1Seq else read2Seq, 
#                 amp$ref1, 
#                 forceQryEnd = QRY_START,
#                 fast = FALSE
#             )
#             swR <- smith_waterman(
#                 unname(rc( if(strand == 0) read2Seq else read1Seq )), 
#                 if(amp$ref2 == "*") amp$ref1 else amp$ref2, 
#                 forceQryEnd = QRY_END,
#                 fast = FALSE
#             )
#             data.table(
#                 molKey = molKey,
#                 dotplotL = list(parseSWDots(swL)),
#                 dotplotR = list(parseSWDots(swR))                
#             )
#         }] 
#     }, mc.cores = env$N_CPU)),
#     by = "molKey"
# )
# #=====================================================================================
