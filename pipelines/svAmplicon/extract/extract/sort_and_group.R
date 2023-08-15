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
#=====================================================================================

#=====================================================================================
message("  loading edges")
edges <- fread(env$INTERIM_FILE, sep = "\t", header = FALSE)
setnames(edges, edgesColumns)
message("  calculating derivative values")
isJunction <- edges[, edgeType != edgeTypes$ALIGNMENT]
isMasked <-edges[, 
    isJunction & 
    (
        edgeType == edgeTypes$MERGE_FAILURE |
        (
            eventSize  < env$MIN_SV_SIZE & 
            insertSize < env$MIN_SV_SIZE 
        )
    )
]
edges[, ":="(
    i = .I,
    molKey = paste(ampliconId, molId, sep = ":"), # used for simplified sorting and grouping
    masked = isMasked,
    node1 = paste(chrom1, strand1, pos1, sep = "/"),
    node2 = paste(chrom2, strand2, pos2, sep = "/")
)]
edges[isJunction, ":="( # junction qualities are the lower of the two flanking alignment qualities
    mapQ     = pmin(mapQ[i - 1],     mapQ[i + 1]),
    alnQ     = pmin(alnQ[i - 1],     alnQ[i + 1]),
    baseQual = pmin(baseQual[i - 1], baseQual[i + 1])
)]
setkey(edges, molKey, i)
#=====================================================================================

#=====================================================================================
message("  counting unique molecules")
nMolecules     <- edges[,                                    length(unique(molKey))]
nMergeFailure  <- edges[edgeType == edgeTypes$MERGE_FAILURE, length(unique(molKey))]
message(paste0("    ", paste(nMolecules,     "unique molecule sequences",     sep = "\t")))
message(paste0("    ", paste(nMergeFailure,  "molecules had merge failures",  sep = "\t")))
#=====================================================================================

#=====================================================================================
message("  collapsing nodes and edges into molecules")
molecules <- edges[, {
    isAlignment <- edgeType == edgeTypes$ALIGNMENT    
    isCandidateJunction <- !isAlignment & !masked
    isMerged <- mergeLevel[1] > 0
    interleave_ <- function(value1, value2){
        paste(c(rbind(value1[isCandidateJunction], value2[isCandidateJunction])), collapse = ":")
    }
    paste_ <- function(value, collapse = ":"){
        paste(value[isCandidateJunction], collapse = collapse)
    }
    min_ <- function(value){
        min(value[isAlignment], na.rm = TRUE)
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

        # junction-level information, concatenated across all candidate, i.e., unmasked junctions
        # for amplicon sequencing, path need not include the outer alignment nodes
        path        = interleave_(node1,    node2),  # node-level info
        chroms      = interleave_(chrom1,   chrom2), # different parsing of the same node-level info as in path
        strands     = interleave_(strand1,  strand2),
        poss        = interleave_(pos1,     pos2),
        edgeTypes   = paste_(edgeType, collapse = ""), # edge-level info
        eventSizes  = paste_(eventSize),
        insertSizes = paste_(insertSize),
        readNs      = paste_(readN1),
        mapQs       = paste_(mapQ),
        alnQs       = paste_(alnQ),
        baseQuals   = paste_(baseQual),

        # read-level information
        minMapQ     = min_(mapQ), # a flag showing just how bad the worst part of a read was
        minAlnQ     = min_(alnQ),
        minBaseQual = min_(baseQual),  
        cigars = paste(cigar[isAlignment], collapse = ":"),   
        seq1  = seq1[1], 
        seq2  = if(!isMerged) seq2[.N]  else "*", # there are only ever one merged or two unmerged reads
        qual1 = qual1[1],
        qual2 = if(!isMerged) qual2[.N] else "*"
    )
}, by = .(molKey)]
#=====================================================================================

#=====================================================================================
message("  aggregating molecules into types, i.e., shared SV paths")
moleculeTypes <- molecules[
    order(
        ampliconId, path, insertSizes, # first columns are for moleculeType grouping (not currently using sequence, just insertSize)
        -isReference, -mergeLevel, -nReadPairs, -minBaseQual, -minMapQ # last columns are   
    )
][, .(
        # molecule-level information
        # ampliconId 
        # path
        # insertSizes
        molKey      = molKey[1],  # one molecule id selected as the index for the type
        molId       = molId[1],        
        isReference = isReference[1],
        nReadPairs  = sum(nReadPairs), # continue summing all original _read pairs_ that has this type 
        nMolecules  = .N, # the number of distinct read pair _sequences_ that had this type  
        mergeLevel  = mergeLevel[1],
        overlap     = overlap[1],
     
        # junction-level information
        chroms      = chroms[1],
        strands     = strands[1],
        poss        = poss[1],
        edgeTypes   = edgeTypes[1],
        eventSizes  = eventSizes[1],   
        readNs      = readNs[1],  
        mapQs       = mapQs[1],
        alnQs       = alnQs[1],
        baseQuals   = baseQuals[1],

        # read-level information
        minMapQ     = max(minMapQ), # thus, take the worst segment score per molecule, but the best of those values across all molecules matching a type
        minAlnQ     = max(minAlnQ),
        minBaseQual = max(minBaseQual),
        cigars      = cigars[1],
        seq1        = seq1[1],
        seq2        = seq2[1],        
        qual1       = qual1[1],
        qual2       = qual2[1]
), by = .(ampliconId, path, insertSizes)][
    order(ampliconId, -nReadPairs) # sort by molecule type frequency
]
#=====================================================================================

#=====================================================================================
message("  tabulating unique SV junctions found in one or more molecules")
molTypeKeys <- moleculeTypes[, molKey]
edges[, isIndexMol := molKey %in% molTypeKeys]
junctions <- edges[edgeType != edgeTypes$ALIGNMENT & masked == FALSE, .( # only tabulate candidate junction, not merge failures our under-sized SVs
    # ampliconId
    # node1
    # node2
    # insertSize
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
    mergeLevel  = max(mergeLevel),
    molTypeKeys = list(unique(molKey[isIndexMol]))
), by = .(ampliconId, node1, node2, insertSize)] # here, include the junction insertion in the definiton, i.e. D5I0 and D5I22 are different junctions, but not the bases themselves
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
