# sort and group molecules based on the SVs they carry

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("  initializing")
library(data.table)
library(parallel)
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
    "cigar1",
    "readN1",
    "chrom2",
    "strand2",
    "pos2",
    "seq2",
    "qual2",    
    "cigar2",
    "readN2",
    "edgeType",
    "mapQ",
    "size",
    "insSize",
    "ampliconId", 
    "mergeLevel", 
    "overlap", 
    "isRef", 
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
    PROPER        = "P",
    FUSED_MERGE_FAILURE = "F",
    REJECTED_INDEL = "R",
    FUSED_MERGE_FAILURE_REJECTED_INDEL = "Q"
)
#=====================================================================================

#=====================================================================================
message("  loading edges")
edges <- fread(env$INTERIM_FILE, sep = "\t", header = FALSE)
setnames(edges, edgesColumns)
message("  calculating derivative values")
edges[, ":="(
    i = .I,
    molKey = paste(ampliconId, molId, sep = ":"), # used for simplified sorting and grouping
    edgeClass = ifelse(                           # molecules come to us sorted by molKey
        edgeType %in% c(
            edgeTypes$FUSED_MERGE_FAILURE, # edgeClass collapses all rejected edge types (F,R,Q) into a single alignment (A)
            edgeTypes$REJECTED_INDEL, 
            edgeTypes$FUSED_MERGE_FAILURE_REJECTED_INDEL
        ), 
        edgeTypes$ALIGNMENT, 
        edgeType
    )    
)]
setkey(edges, molKey, i)
#=====================================================================================

#=====================================================================================
message("  counting unique molecules")
nMolecules <- edges[, length(unique(molKey))]
nMergeFailure  <- edges[edgeType %in% c(edgeTypes$FUSED_MERGE_FAILURE, edgeTypes$FUSED_MERGE_FAILURE_REJECTED_INDEL), length(unique(molKey))]
nRejectedIndel <- edges[edgeType %in% c(edgeTypes$REJECTED_INDEL,      edgeTypes$FUSED_MERGE_FAILURE_REJECTED_INDEL), length(unique(molKey))]
nBoth <- edges[edgeType == edgeTypes$FUSED_MERGE_FAILURE_REJECTED_INDEL, length(unique(molKey))]
message(paste0("    ", paste(nMolecules,     "unique molecule sequences",     sep = "\t")))
message(paste0("    ", paste(nMergeFailure,  "molecules had merge failures",  sep = "\t")))
message(paste0("    ", paste(nRejectedIndel, "molecules had rejected indels", sep = "\t")))
message(paste0("    ", paste(nBoth,          "molecules had both merge failures and rejected indels", sep = "\t")))
#=====================================================================================

#=====================================================================================
message("  collapsing nodes and edges into molecules")
edges[, ":="(
    nodePair = paste( # create a unique identifier for an edge, NOT including insSize yet, only for path tracing
        paste(chrom1, strand1, pos1, sep = "/"),
        paste(chrom2, strand2, pos2, sep = "/"),
        sep = ":"          
    )
)]
molecules <- edges[, {
    isAlignment <- edgeClass == edgeTypes$ALIGNMENT
    isSameRead <- readN1 == readN2
    isRead1 <- readN1 == 1
    isMerged <- mergeLevel[1] == 0
    aligmentInfo <- function(value1, value2){
        paste0(na.omit(ifelse(
            isAlignment, 
            ifelse(
                isSameRead, 
                ifelse(isRead1, value1, value2), # read 2 reports end positions, for consistency with merge failures
                paste(value1, value2, sep = ":") # properly handle merge failures, i.e., 1:2 fusions, by reporting both alignments
            ),
            NA
        )), collapse = ":")
    }
    junctionInfo <- function(value, collapse = ":"){
        paste0(na.omit(ifelse(
            isAlignment, 
            NA,
            value
        )), collapse = collapse)
    }
    .(
        # molecule-level information, the same for every edge
        ampliconId = ampliconId[1], 
        molId = molId[1],  
        isRef =  as.logical(isRef[1]),    
        nReadPairs = nReadPairs[1],    
        mergeLevel = mergeLevel[1], 
        overlap = overlap[1],

        # alignment-level information
        chroms = aligmentInfo(chrom1,  chrom2),
        strands= aligmentInfo(strand1, strand2),
        poss   = aligmentInfo(pos1,    pos2),
        cigars = aligmentInfo(cigar1,  cigar2),
        readNs = aligmentInfo(readN1,  readN2),
        mapQs       = aligmentInfo(mapQ,     mapQ),
        baseQuals   = aligmentInfo(baseQual, baseQual),
        minMapQ     = min(mapQ,     na.rm = TRUE),
        minBaseQual = min(baseQual, na.rm = TRUE),

        # junction-level information
        pathClass = junctionInfo(edgeClass, collapse = ""),
        path      = junctionInfo(nodePair),
        sizes     = junctionInfo(size),
        insSizes  = junctionInfo(insSize),
        nodePairs = list(as.character(na.omit(ifelse( # a concatenation of internal nodePairs only; zero-length string when no SV found
            edgeClass == edgeTypes$ALIGNMENT, 
            NA, 
            nodePair
        )))),

        # read-level information
        seq1  = seq1[1], 
        seq2  = if(isMerged) seq2[.N]  else "*", # there are only ever one merged or two unmerged reads
        qual1 = qual1[1],
        qual2 = if(isMerged) qual2[.N] else "*"
    )
}, by = .(molKey)]
#=====================================================================================

#=====================================================================================
message("  aggregating molecules into types, i.e., shared SV paths")
moleculeTypes <- molecules[
    order(
        ampliconId, path, insSizes,      # first columns are for next grouping
        -isRef, -mergeLevel, -nReadPairs # last columns are for selecting the index molecule   
    )
][,
    .(
        # molecule-level information
        # ampliconId 
        # path
        # insSizes
        molKey      = molKey[1],
        molId       = molId[1],  # one molecule id selected as the index for the type        
        isRef       = isRef[1],
        nReadPairs  = sum(nReadPairs), # continue summing all original read pairs that has this type 
        nMols       = .N, # the number of distinct read pair _sequences_ that had this type        
        mergeLevel  = mergeLevel[1],
        overlap     = overlap[1],
     
        # alignment-level information
        chroms      = chroms[1],
        strands     = strands[1],
        poss        = poss[1],
        cigars      = cigars[1],
        readNs      = readNs[1],        
        mapQs       = mapQs[1],
        baseQuals   = baseQuals[1],
        minMapQ     = max(minMapQ), # thus, take the worst segment score per molecule, but the best of those values across all molecules matching a type
        minBaseQual = max(minBaseQual),

        # junction-level information
        pathClass   = pathClass[1],
        nodePairs   = nodePairs[1],
        sizes       = sizes[1],        

        # read-level information
        seq1        = seq1[1],
        seq2        = seq2[1],        
        qual1       = qual1[1],
        qual2       = qual2[1]
    ), by = .(ampliconId, path, insSizes) # group
][
    order(ampliconId, -nReadPairs) # sort by molecule type frequency
]
#=====================================================================================

#=====================================================================================
message("  tabulating unique SV junctions found in one or more molecules")
molTypeKeys <- moleculeTypes[, molKey]
edges[, isIndexMol := molKey %in% molTypeKeys]
junctions <- edges[edgeClass != edgeTypes$ALIGNMENT, .(
    # ampliconId
    # nodePair
    # insSize
    nReadPairs = sum(nReadPairs),
    nMol = length(unique(molKey)), # beware of possibility of rare molecules with a recurring junction
    nMolTypes = length(unique(molKey[isIndexMol])), 
    edgeClass = edgeClass[1],
    mergeLevel = max(mergeLevel),
    size = size[1],
    mapQ = max(mapQ),
    molTypeKeys = list(unique(molKey[isIndexMol]))
), by = .(ampliconId, nodePair, insSize)] # here, include the junctional insertion in the definiton, i.e. D5I0 and D5I22 are different junctions
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
