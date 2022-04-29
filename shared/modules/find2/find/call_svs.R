
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
        'GENOME_CHROMS',
        'GENOME_FASTA',
        'CHROM_FASTA_DIR'        
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
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) ########################
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general', 'IUPAC', 'smith_waterman'))
sourceScripts(file.path(rUtilDir, 'genome'),     c('general', 'chroms', 'targets', 'faidx'))
sourceScripts(file.path(env$ACTION_DIR, 'find'), c(
    'column_definitions', 'parse_groups', 'purge_duplicates', 
    'assign_references', 'get_outer_clips', 'analyze_junctions'
))
#-------------------------------------------------------------------------------------
# initialize samples
message("initializing samples")
env$SAMPLES <- strsplit(env$SAMPLES, "\\s+")[[1]]
SAMPLES <- as.list(seq_along(env$SAMPLES)) # create sample indices
names(SAMPLES) <- env$SAMPLES
env$MAX_TLENS <- as.integer(strsplit(env$MAX_TLENS, "\\s+")[[1]])
MAX_MAX_TLEN <- max(env$MAX_TLENS)
# MAX_TLENS <- as.list(env$MAX_TLENS)
# names(MAX_TLENS) <- env$SAMPLES
#=====================================================================================

#=====================================================================================
# load and preprocess the called and grouped SVs across all samples
#-------------------------------------------------------------------------------------
message("loading and parsing molecules")
jxnMols <- as.data.table(read.table(
    'stdin',
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = names(find$working),
    colClasses = unname(unlist(find$working))
))
jxnMols[, c('chrom1', 'side1', 'pos1')] <- unpackNodeNames(jxnMols$NODE_1)
jxnMols[, c('chrom2', 'side2', 'pos2')] <- unpackNodeNames(jxnMols$NODE_2)
jxnMols[, ':='(
    jxnName       = paste(NODE_1, NODE_2, sep = ","), # the junction edge called by a molecule (could be >1 per molecule) # nolint
    jxnKey        = paste(SAMPLES[SAMPLE], MOL_ID, JXN_N, sep = ":"), # sample-level ID for the source junction edge (each with 2 nodes) # nolint
    svIndex       = '', # SV identifier shared between all samples
    sampleSvIndex = '', # SV identifier unique to one sample
    AMBIGUOUS     = 0L, # AMBIGUOUS gap molecules are consistent with more than one split and present in the table more than once # nolint
    DOWNSAMPLED   = 0L, # DOWNSAMPLED molecules were originally part of a larger molecule set, some of which were dropped # nolint
    N_COLLAPSED   = 1L, # the number of input molecules collapsed as duplicates into this remaining molecule
    IS_REFERENCE  = 0L  # reference molecules acted as the primary guide to characterizing an SV
)]
#=====================================================================================

#=====================================================================================
# apply initial molecule-level filters
#-------------------------------------------------------------------------------------
message("applying SV size filter to junction molecules")
if(env$SV_SIZE_FACTOR > 0) env$MIN_SV_SIZE <- env$SV_SIZE_FACTOR * MAX_MAX_TLEN
if(env$MIN_SV_SIZE > 0) jxnMols <- jxnMols[JXN_TYPE == "T" | abs(pos2 - pos1) >= env$MIN_SV_SIZE]
#=====================================================================================

#=====================================================================================
# break molecule continuity groups into SV calls
#-------------------------------------------------------------------------------------
message("parsing junction continuity groups into SVs")
jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(groupIndex)], parseContinuityGroup, mc.cores = env$N_CPU))

message("marking ambiguous source molecules")
setkey(jxnMols, jxnKey)
ambiguousJxnKeys <- jxnMols[, .N, by = jxnKey][N > 1, jxnKey]
jxnMols[ambiguousJxnKeys, AMBIGUOUS := 1L]
#=====================================================================================

#=====================================================================================
# purge duplicates from within the molecules for each SV in each sample
#-------------------------------------------------------------------------------------
message("aggregating source molecule duplicates")
setkey(jxnMols, sampleSvIndex)
jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(sampleSvIndex)], purgeDuplicateMolecules, mc.cores = env$N_CPU))
#=====================================================================================

#=====================================================================================
# assign one molecule as the reference molecule for further characterization of each SV junction
#-------------------------------------------------------------------------------------
message("assigning SV reference molecules")
setkey(jxnMols, svIndex)
jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], assignReferenceMolecule, mc.cores = env$N_CPU))
#=====================================================================================

#=====================================================================================
# add reliable outer clip nodes as SV evidence and use to knit together gap junctions
#-------------------------------------------------------------------------------------
message("collecting candidate matching outer clipped nodes")
refNodesFile <- printReferenceNodes()
setkey(jxnMols, svIndex) # key does not persist after do.call(rbind, mclapply...
flipGuidance <- jxnMols[IS_REFERENCE == 1, .(
    svIndex     = svIndex,
    isInversion = side1 == side2, 
    flipNode    = ifelse(side1 == 'R', 1L, 2L)
)]
jxnMols <- rbind(jxnMols, do.call(rbind, getOuterClipEvidence())) # one sample at a time with internal parallel actions 
unlink(refNodesFile)

message("attempting to merge gap junctions using clipped nodes")
jxnMols[IS_REFERENCE == 1, ':='(
    JXN_SEQ = ifelse(NODE_CLASS == nodeClasses$SPLIT, SEQ_1, "*"),
    MERGE_LEN = 0L
)]
setkey(jxnMols, svIndex)
jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], mergeGapJunctions, mc.cores = env$N_CPU))

message("purging invalid/untrustworthy clip evidence")
setkey(jxnMols, svIndex)
jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], purgeInvalidClips, mc.cores = env$N_CPU))
#=====================================================================================

#=====================================================================================
# characterize the final set of SV junctions and prepare for app
#-------------------------------------------------------------------------------------
# initialize genome
message("initializing genome sequence retrieval")
setCanonicalChroms()
write(
    paste0('CHROMS: ', paste(canonicalChroms, collapse = " ")), 
    file = paste(env$FIND_PREFIX, "metadata", "yml", sep = "."),
    append = TRUE
)
loadFaidx() 
faidx_padding <- as.integer(MAX_MAX_TLEN * 1.2) # sufficient to contain any source molecule

# initialize and store parsed target regions (if any)
loadTargetRegions()
write.table(
    if(!is.null(targetRegions)) targetRegions$bed else "NA", 
    paste(env$FIND_PREFIX, "target_regions", "bed", sep = "."), 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)

message("characterizing final SV junction calls")
setkey(jxnMols, svIndex)
svCalls <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], characterizeSvJunction, mc.cores = env$N_CPU))

message("counting SVs by sample")
svCalls <- merge(svCalls, dcast(jxnMols, svIndex ~ SAMPLE, length), by.x = 'SV_ID', by.y = "svIndex")

message("setting SV molecule plot positions")
jxnMols[, TARGET_POS_1 := mcmapply(getTargetRegionI, chrom1, OUT_POS1, mc.cores = env$N_CPU)]
jxnMols[, TARGET_POS_2 := mcmapply(getTargetRegionI, chrom2, OUT_POS2, mc.cores = env$N_CPU)]
#=====================================================================================

#=====================================================================================
# save results and print some summary statistics
#-------------------------------------------------------------------------------------
message("writing SV summary table")
outFile <- paste(env$FIND_PREFIX, 'structural_variants', 'gz', sep = ".")
fwrite(
    svCalls, 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = TRUE, 
    compress = "gzip"
)
outFile <- paste(env$FIND_PREFIX, 'structural_variants', 'rds', sep = ".")
saveRDS(
    svCalls, 
    file = outFile
)

message("writing molecule evidence table")
outFile <- paste(env$FIND_PREFIX, 'junction_molecules', 'gz', sep = ".")
jxnMols[, SV_ID := svIndex]
jxnMols <- jxnMols[, .SD, .SDcols = names(find$junction_molecules)]
fwrite(
    jxnMols, 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = TRUE, 
    compress = "gzip"
)
outFile <- paste(env$FIND_PREFIX, 'junction_molecules', 'rds', sep = ".")
saveRDS(
    jxnMols, 
    file = outFile
)
#=====================================================================================
