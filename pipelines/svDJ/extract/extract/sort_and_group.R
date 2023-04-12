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
        'JUNCTIONS_FILE'
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
#-------------------------------------------------------------------------------------
# source R scripts
# sourceScripts(rUtilDir, 'utilities')
# rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
# sourceScripts(file.path(rUtilDir, 'sequence'), c('general','IUPAC','smith_waterman'))
#-------------------------------------------------------------------------------------
# set some constants
edgesColumns <- c(
    "molId", # will be become the type id once an index is selected
    "chrom1", 
    "strand1",
    "pos1",
    "seq1",
    "qual1",    
    "cigar1",
    "chrom2",
    "strand2",
    "pos2",
    "seq2",
    "qual2",    
    "cigar2",
    "edgeType",
    "mapQ",
    "size",
    "insSize",
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
    PROPER        = "P",
)
#=====================================================================================

#=====================================================================================
message("  loading edges")
edges <- fread(env$INTERIM_FILE, sep = "\t", header = FALSE)
setnames(edges, edgesColumns)
message("  calculating derivative values")
edges[, ":="(i = .I)]
setkey(edges, molId, i)
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
    isAlignment <- edgeType == edgeTypes$ALIGNMENT
    aligmentInfo <- function(value){
        paste0(na.omit(ifelse(isAlignment, value, NA)), collapse = ":")
    }
    junctionInfo <- function(value, collapse = ":"){
        paste0(na.omit(ifelse(isAlignment, NA, value)), collapse = collapse)
    }
    .(
        # alignment-level information
        chroms = aligmentInfo(chrom1),
        strands= aligmentInfo(strand1),
        poss   = aligmentInfo(pos1),
        cigars = aligmentInfo(cigar1),
        mapQs       = aligmentInfo(mapQ),
        baseQuals   = aligmentInfo(baseQual),
        minMapQ     = min(mapQ,     na.rm = TRUE),
        minBaseQual = min(baseQual, na.rm = TRUE),

        # junction-level information
        pathClass = junctionInfo(edgeType, collapse = ""),
        path      = junctionInfo(nodePair),
        sizes     = junctionInfo(size),
        insSizes  = junctionInfo(insSize),
        nodePairs = list(as.character(na.omit(ifelse(edgeType == edgeTypes$ALIGNMENT, NA, nodePair)))),# a concatenation of internal nodePairs only; zero-length string when no SV found

        # read-level information
        seq  = seq1[1], 
        qual = qual1[1],
    )
}, by = .(molId)]
#=====================================================================================

#=====================================================================================
message("  aggregating molecules into types, i.e., shared junction paths")
moleculeTypes <- molecules[
    order(
        path, insSizes,        # first columns are for next grouping
        -minBaseQual, -minMapQ # last columns are for selecting the index molecule   
    )
][,
    .(
        # molecule-level information
        # path
        # insSizes
        molId       = molId[1],  # one molecule id selected as the index for the type        
        nMols       = .N, # the number of distinct molecules that had this type        
     
        # alignment-level information
        chroms      = chroms[1],
        strands     = strands[1],
        poss        = poss[1],
        cigars      = cigars[1],      
        mapQs       = mapQs[1],
        baseQuals   = baseQuals[1],
        minMapQ     = max(minMapQ), # thus, take the worst segment score per molecule, but the best of those values across all molecules matching a type
        minBaseQual = max(minBaseQual),

        # junction-level information
        pathClass   = pathClass[1],
        nodePairs   = nodePairs[1],
        sizes       = sizes[1],        

        # read-level information
        seq        = seq[1],     
        qual       = qual[1]
    ), by = .(path, insSizes) # group
][
    order(-nMols) # sort by molecule type frequency
]
#=====================================================================================

#=====================================================================================
message("  tabulating unique junctions found in one or more molecules")
molTypeIds <- moleculeTypes[, molId]
edges[, isIndexMol := molId %in% molTypeIds]
junctions <- edges[edgeType != edgeTypes$ALIGNMENT, .(
    # nodePair
    # insSize
    nMol = length(unique(molId)), # beware of possibility of rare molecules with a recurring junction
    nMolTypes = length(unique(molId[isIndexMol])), 
    edgeType = edgeType[1],
    size = size[1],
    mapQ = max(mapQ),
    molTypeIds = list(unique(molId[isIndexMol]))
), by = .(nodePair, insSize)] # here, include the junctional insertion in the definiton, i.e. D5I0 and D5I22 are different junctions
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
