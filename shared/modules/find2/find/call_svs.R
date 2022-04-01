
# compare structural variants between samples and annotate SV calls with any matching samples

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(parallel)
library(data.table)
library(jsonlite)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'DATA_NAME',
        'ACTION_DIR',
        'MODULES_DIR',
        'COMPILE_PREFIX',
        'FIND_PREFIX',
        'FIND_MODE',
        'SHM_DIR_WRK',
        'SAMPLES',
        'MAX_TLENS',
        'GENOME',
        'GENOME_CHROMS'
    ),
    integer = c(
        'N_CPU',

        'MIN_MAPQ_ONE',
        'MIN_MAPQ_BOTH',

        'MIN_SV_SIZE',
        'SV_SIZE_FACTOR',

        'PURGE_DISTANCE',
        'PURGE_LIMIT',
        'MIN_COVERAGE',
        'MIN_MERGE_OVERLAP'
    ),
    double = c(
        'MIN_MERGE_DENSITY'
    )   
))
# writeEnvJson(env$FIND_PREFIX) # for use by Shiny app
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general', 'IUPAC', 'smith_waterman'))
sourceScripts(file.path(rUtilDir, 'genome'),     c('general', 'chroms', 'faidx'))
sourceScripts(file.path(env$ACTION_DIR, 'find'), c(
    'column_definitions', 'node_retrieval', 
    'parse_groups', 'purge_duplicates', 'assign_references', 'get_outer_clips', 'analyze_junctions'
))
# , 'network'
#-------------------------------------------------------------------------------------
# initialize samples
env$SAMPLES <- strsplit(env$SAMPLES, "\\s+")[[1]]
SAMPLES <- as.list(1:length(env$SAMPLES)) # create sample indices
names(SAMPLES) <- env$SAMPLES
env$MAX_TLENS <- as.integer(strsplit(env$MAX_TLENS,"\\s+")[[1]])
MAX_MAX_TLEN <- max(env$MAX_TLENS)
MAX_TLENS <- as.list(env$MAX_TLENS)
names(MAX_TLENS) <- env$SAMPLES
faidx_padding <- round(MAX_MAX_TLEN * 1.2, 0) # sufficient to contain any source molecule span
#-------------------------------------------------------------------------------------
# initialize genome
setCanonicalChroms()
# loadFaidx(env$SHM_DIR_WRK)
#=====================================================================================

#=====================================================================================
# load and preprocess the called and grouped SVs across all samples
#-------------------------------------------------------------------------------------
jxnMols <- as.data.table(read.table(
    'stdin',
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = names(compile$working),
    colClasses = unname(unlist(compile$working))
))
#-------------------------------------------------------------------------------------
jxnMols[, ':='(
    AMBIGUOUS = 0, # ambiguous gap molecules are consistent with more than one split
    N_COLLAPSED = 1, # the number of input molecules collapsed as duplicates into the remaining molecule
    jxnName   = paste(NODE_1, NODE_2, sep = ","), # the junction edge called by a molecule (could be >1 per molecule)
    jxnKey    = paste(SAMPLES[SAMPLE], MOL_ID, JXN_N, sep = ":") # ID for the source junction edge (each with 2 nodes)
)]
#-------------------------------------------------------------------------------------
jxnMols[, c('chrom1', 'side1', 'pos1')] <- unpackNodeNames(jxnMols$NODE_1)
jxnMols[, c('chrom2', 'side2', 'pos2')] <- unpackNodeNames(jxnMols$NODE_2)
if(env$SV_SIZE_FACTOR > 0) env$MIN_SV_SIZE <- env$SV_SIZE_FACTOR * MAX_MAX_TLEN
if(env$MIN_SV_SIZE > 0) jxnMols <- jxnMols[JXN_TYPE == "T" | abs(pos2 - pos1) >= env$MIN_SV_SIZE]
#=====================================================================================

# #=====================================================================================
# # continue breaking molecule continuity groups into SV calls
# #-------------------------------------------------------------------------------------
# message(nrow(jxnMols))
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(groupIndex)], parseContinuityGroup, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(groupIndex)], parseContinuityGroup))
# ambiguousJxnKeys <- jxnMols[, .N, by = jxnKey][N > 1, jxnKey]
# jxnMols[jxnKey %in% ambiguousJxnKeys, AMBIGUOUS := 1] # mark ambiguous gap molecules present more than once in table
# message(nrow(jxnMols))
# message(jxnMols[, length(unique(sampleSvIndex))])
# #=====================================================================================

# #=====================================================================================
# # purge duplicates from within the molecules for each SV in each sample
# #-------------------------------------------------------------------------------------
# message("purging molecule duplicates")
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(sampleSvIndex)], purgeDuplicateMolecules, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(sampleSvIndex)], purgeDuplicateMolecules))
# message(nrow(jxnMols))
# message(jxnMols[, length(unique(svIndex))])
# #=====================================================================================

##############################
tmpFile <- paste(env$FIND_PREFIX, "DEVELOP.RDS", sep=".")
# saveRDS(jxnMols, file = tmpFile)
# jxnMols <- readRDS(tmpFile)

# #=====================================================================================
# # assign one molecule as the reference molecule for further characterization of each SV junction
# #-------------------------------------------------------------------------------------
# message("assigning SV reference molecules")
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], assignReferenceMolecule, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(svIndex)[1:1000]], assignReferenceMolecule))

# ##############################
# # saveRDS(jxnMols, file = tmpFile)
# jxnMols <- readRDS(tmpFile)

# refMols <- jxnMols[, IS_REFERENCE == 1]
# refNodes <- printReferenceNodes()
# # #=====================================================================================

# #=====================================================================================
# # add out clip nodes as SV evidence, one sample at a time with internal parallel actions
# #-------------------------------------------------------------------------------------
# message("adding outer clip nodes")
# flipGuidance <- jxnMols[refMols, .(
#     svIndex     = svIndex,
#     isInversion = side1 == side2, # correct for svCaptur, is it correct for svWGS
#     flipNode    = ifelse(side1 == 'R', 1L, 2L)
# )]
# setkey(flipGuidance, svIndex)
# jxnMols <- rbind(jxnMols, do.call(rbind, getOuterClipEvidence()))
# #=====================================================================================

tmpFile <- paste(env$FIND_PREFIX, "DEVELOP.WITH_CLIPS.RDS", sep=".")
# saveRDS(jxnMols, file = tmpFile)
jxnMols <- readRDS(tmpFile)

#=====================================================================================
# assemble and characterize final junction calls
#-------------------------------------------------------------------------------------
message("characterizing SV junctions")
jxnMols[IS_REFERENCE == 1, ':='(
    JXN_SEQ = ifelse(NODE_CLASS == nodeClasses$SPLIT, SEQ_1, NA),
    MERGE_LEN = 0,
    RECONSTRUCTED = 0
)]
# setkey(jxnMols, svIndex)
jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], mergeGapJunctions, mc.cores = env$N_CPU))
# jxnMols <- do.call(rbind, lapply(jxnMols[, unique(svIndex)[1:1000]], mergeGapJunctions))

str(jxnMols[RECONSTRUCTED != 0])
stop("ccc  XXXXXXXXXXXXXXXXX")
#=====================================================================================




#=====================================================================================
# save results and print some summary statistics
#-------------------------------------------------------------------------------------
message("writing SV summary table")
outFile <- paste(env$COMPARE_PREFIX, 'structural_variants', 'gz', sep = ".")
outFile <- gzfile(outFile, "w")
write.table(
    svs[, .SD, .SDcols = names(compare$structural_variants)], 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = FALSE
)
close(outFile)
#-------------------------------------------------------------------------------------
reportStat(nrow(svs),                "SVs called across all samples")
reportStat(nrow(svs[N_MATCHES > 0]), "SVs matched at least one other sample")
#=====================================================================================
