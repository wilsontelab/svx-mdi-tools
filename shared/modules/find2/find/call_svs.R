
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
sourceScripts(file.path(rUtilDir, 'genome'),     c('general', 'chroms', 'targets', 'faidx'))
sourceScripts(file.path(env$ACTION_DIR, 'find'), c(
    'column_definitions', 'parse_groups', 'purge_duplicates', 
    'assign_references', 'get_outer_clips', 'analyze_junctions'
))
#-------------------------------------------------------------------------------------
# initialize samples
message("initializing samples")
env$SAMPLES <- strsplit(env$SAMPLES, "\\s+")[[1]]
SAMPLES <- as.list(1:length(env$SAMPLES)) # create sample indices
names(SAMPLES) <- env$SAMPLES
env$MAX_TLENS <- as.integer(strsplit(env$MAX_TLENS,"\\s+")[[1]])
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
    col.names = names(compile$working),
    colClasses = unname(unlist(compile$working))
))
jxnMols[, c('chrom1', 'side1', 'pos1')] <- unpackNodeNames(jxnMols$NODE_1)
jxnMols[, c('chrom2', 'side2', 'pos2')] <- unpackNodeNames(jxnMols$NODE_2)
jxnMols[, ':='(
    jxnName       = paste(NODE_1, NODE_2, sep = ","), # the junction edge called by a molecule (could be >1 per molecule)
    jxnKey        = paste(SAMPLES[SAMPLE], MOL_ID, JXN_N, sep = ":"), # sample-level ID for the source junction edge (each with 2 nodes)
    svIndex       = '', # SV identifier shared between all samples
    sampleSvIndex = '', # SV identifier unique to one sample
    AMBIGUOUS     = 0L, # AMBIGUOUS gap molecules are consistent with more than one split and present in the table more than once
    DOWNSAMPLED   = 0L, # DOWNSAMPLED molecules were originally part of a larger molecule set, some of which were dropped
    N_COLLAPSED   = 1L, # the number of input molecules collapsed as duplicates into this remaining molecule
    IS_REFERENCE  = 0L  # reference molecules acted as the primary guide to characterizing an SV
)]
#=====================================================================================

# #=====================================================================================
# # apply initial molecule-level filters
# #-------------------------------------------------------------------------------------
# message("applying SV size filter to junction molecules")
# if(env$SV_SIZE_FACTOR > 0) env$MIN_SV_SIZE <- env$SV_SIZE_FACTOR * MAX_MAX_TLEN
# if(env$MIN_SV_SIZE > 0) jxnMols <- jxnMols[JXN_TYPE == "T" | abs(pos2 - pos1) >= env$MIN_SV_SIZE]
# #=====================================================================================

# #=====================================================================================
# # continue breaking molecule continuity groups into SV calls
# #-------------------------------------------------------------------------------------
# message("parsing continuity groups into SVs")
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(groupIndex)], parseContinuityGroup, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(groupIndex)], parseContinuityGroup))

# message("marking ambiguous molecules")
# setkey(jxnMols, jxnKey)
# ambiguousJxnKeys <- jxnMols[, .N, by = jxnKey][N > 1, jxnKey]
# jxnMols[ambiguousJxnKeys, AMBIGUOUS := 1]
# #=====================================================================================

# #=====================================================================================
# # purge duplicates from within the molecules for each SV in each sample
# #-------------------------------------------------------------------------------------
# message("purging molecule duplicates")
# setkey(jxnMols, sampleSvIndex)
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(sampleSvIndex)], purgeDuplicateMolecules, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(sampleSvIndex)], purgeDuplicateMolecules))
# #=====================================================================================

# ##############################
# tmpFile <- paste(env$FIND_PREFIX, "DEVELOP_1.RDS", sep=".")
# saveRDS(jxnMols, file = tmpFile)
# jxnMols <- readRDS(tmpFile)

# #=====================================================================================
# # assign one molecule as the reference molecule for further characterization of each SV junction
# #-------------------------------------------------------------------------------------
# message("assigning SV reference molecules")
# setkey(jxnMols, svIndex)
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], assignReferenceMolecule, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(svIndex)[1:1000]], assignReferenceMolecule))
# #=====================================================================================

# ##############################
# tmpFile <- paste(env$FIND_PREFIX, "DEVELOP_2.RDS", sep=".")
# saveRDS(jxnMols, file = tmpFile)
# jxnMols <- readRDS(tmpFile)

# #=====================================================================================
# # add out clip nodes as SV evidence, one sample at a time with internal parallel actions
# #-------------------------------------------------------------------------------------
# message("collecting candidate matching outer clip nodes")
# refNodesFile <- printReferenceNodes()
# setkey(jxnMols, svIndex) # key does not persist after do.call(rbind, mclapply...
# flipGuidance <- jxnMols[IS_REFERENCE == 1, .(
#     svIndex     = svIndex,
#     isInversion = side1 == side2, 
#     flipNode    = ifelse(side1 == 'R', 1L, 2L)
# )]
# jxnMols <- rbind(jxnMols, do.call(rbind, getOuterClipEvidence()))
# unlink(refNodesFile)
# #=====================================================================================

# ##############################
# tmpFile <- paste(env$FIND_PREFIX, "DEVELOP_3.RDS", sep=".")
# saveRDS(jxnMols, file = tmpFile)
# jxnMols <- readRDS(tmpFile)

# #=====================================================================================
# # assemble and characterize final junction calls
# #-------------------------------------------------------------------------------------
# message("attempting to merge gap junctions")
# jxnMols[IS_REFERENCE == 1, ':='(
#     JXN_SEQ = ifelse(NODE_CLASS == nodeClasses$SPLIT, SEQ_1, "*"),
#     MERGE_LEN = 0L
# )]
# setkey(jxnMols, svIndex)
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], mergeGapJunctions, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(svIndex)[1:1000]], mergeGapJunctions))

# message("purging invalid/untrustworthy clip evidence")
# setkey(jxnMols, svIndex)
# jxnMols <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], purgeInvalidClips, mc.cores = env$N_CPU))
# # jxnMols <- do.call(rbind, lapply(jxnMols[, unique(svIndex)[1:1000]], purgeInvalidClips))

# ##############################
# tmpFile <- paste(env$FIND_PREFIX, "DEVELOP_4.RDS", sep=".")
# # saveRDS(jxnMols, file = tmpFile)
# jxnMols <- readRDS(tmpFile)

# initialize genome
message("initializing genome sequence retrieval")
setCanonicalChroms()
loadFaidx(env$SHM_DIR_WRK)
faidx_padding <- round(MAX_MAX_TLEN * 1.2, 0) # sufficient to contain any source molecule span

# initialize target regions (if any)
loadTargetRegions()
write.table(
    if(!is.null(targetRegions)) targetRegions$bed else "NA", 
    paste(env$FIND_PREFIX, "target_regions", "bed", sep="."), 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)

# message("characterizing final SV junction calls")
# setkey(jxnMols, svIndex)
# svCalls <- do.call(rbind, mclapply(jxnMols[, unique(svIndex)], characterizeSvJunction, mc.cores = env$N_CPU))
# # svCalls <- do.call(rbind, lapply(jxnMols[, unique(svIndex)], characterizeSvJunction))

# message("counting SVs by sample")
# svCalls <- merge(svCalls, dcast(jxnMols, svIndex ~ SAMPLE, length), by.x = 'SV_ID', by.y = "svIndex")

##############################
jxnMolsFile <- paste(env$FIND_PREFIX, "DEVELOP_5_jxnMols.RDS", sep=".")
svCallsFile <- paste(env$FIND_PREFIX, "DEVELOP_5_svCalls.RDS", sep=".")
# saveRDS(jxnMols, file = jxnMolsFile)
# saveRDS(svCalls, file = svCallsFile)
jxnMols <- readRDS(jxnMolsFile)
svCalls <- readRDS(svCallsFile)

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

message("writing molecule evidence table")
outFile <- paste(env$FIND_PREFIX, 'junction_molecules', 'gz', sep = ".")
fwrite(
    jxnMols, # TODO, restrict these output columns jxnMols[, .SD, .SDcols = names(compare$junction_molecules)], 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = TRUE, 
    compress = "gzip"
)
#=====================================================================================

#  $ NODE_1        : chr  "16:R:78523150" "16:R:78523183" "16:R:78523199" "16:R:78523164" ...
#  $ CLIP_LEN_1    : int  149 0 0 0 109 31 31 31 0 0 ...
#  $ CLIP_SEQ_1    : chr  "GATCATGCCACTGCACTCCAGCTTGGGCAAAAGAGCAAGACTCATCTCAAACAAACAAATAAATAAATAAATAAATAAATAAATAAATGAAATCTCTAGTCAACTAAACTT"| __truncated__ "*" "*" "*" ...
#  $ FLAG_1        : int  2128 81 81 81 80 81 81 81 81 81 ...
#  $ POS_1         : int  78523150 78523183 78523199 78523164 78523150 78523150 78523150 78523150 78523216 78523183 ...
#  $ MAPQ_1        : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ CIGAR_1       : chr  "120M149S" "151M" "151M" "151M" ...
#  $ SEQ_1         : chr  "CACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACGTGCCTTTAAATGAAGATGGTAATTATTACTATTACTGCCACCACCTTTGATAAATTTAGGTCTG"| __truncated__ "GCAGGGCATCCAGTTTACAGCTGTCTACCTTACCTTGAATTTATATATACTGGCACATATTAAACACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACG"| __truncated__ "GCATTTGCTGGGTACAGCAGGGCATCCAGTTTACAGCTGTCTACCTTACCTTGAATTTATATATACTGGCACATATTAAACACTATGAAAAATTAATATTTTATTTAGGTC"| __truncated__ "GCTGTCTACCTTACCTTGAATTTATATATACTGGCACATATTAAACACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACGTGCCTTTAAATGAAGATGG"| __truncated__ ...
#  $ ALN_N_1       : int  1 1 1 1 2 2 2 2 1 1 ...
#  $ UMI_1         : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ NODE_CLASS    : int  1 0 0 0 1 1 1 1 0 0 ...
#  $ JXN_TYPE      : chr  "I" "I" "I" "I" ...
#  $ JXN_N         : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ MOL_ID        : int  455840 455841 455842 714401 714776 714777 714778 714779 714780 518766 ...
#  $ IS_MERGED     : int  1 0 0 0 1 0 0 0 0 0 ...
#  $ IS_DUPLEX     : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ STRAND_COUNT1 : int  10 8 13 5 0 0 0 0 0 14 ...
#  $ STRAND_COUNT2 : int  0 0 0 0 4 2 11 8 2 0 ...
#  $ MOL_CLASS     : chr  "V" "V" "V" "V" ...
#  $ MOL_STRAND    : int  0 0 0 0 1 1 1 1 1 0 ...
#  $ IS_OUTER_CLIP1: int  0 0 0 0 0 0 0 0 0 0 ...
#  $ IS_OUTER_CLIP2: int  0 0 0 0 0 0 0 0 0 0 ...
#  $ TARGET_CLASS  : chr  "TT" "TT" "TT" "TT" ...
#  $ SHARED_PROPER : int  2 2 2 2 2 2 2 2 2 2 ...
#  $ OUT_POS1      : int  78523269 78523333 78523349 78523314 78523269 78523269 78523269 78523269 78523366 78523333 ...
#  $ OUT_POS2      : int  78524487 78524494 78524487 78524535 78524447 78524487 78524494 78524549 78524691 78524487 ...
#  $ SAMPLE        : chr  "RO3306_0.2APH_Colch_2APH_M" "RO3306_0.2APH_Colch_2APH_M" "RO3306_0.2APH_Colch_2APH_M" "RO3306_0.2APH_Colch_M" ...
#  $ NODE_2        : chr  "16:R:78524315" "16:R:78524332" "16:R:78524325" "16:R:78524373" ...
#  $ CLIP_LEN_2    : int  108 0 0 0 108 108 108 108 0 0 ...
#  $ CLIP_SEQ_2    : chr  "CACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACGTGCCTTTAAATGAAGATGGTAATTATTACTATTACTGCCACCACCTTTGATAAATTTAGGT" "*" "*" "*" ...
#  $ FLAG_2        : int  128 161 161 161 2176 2177 2177 2177 161 161 ...
#  $ POS_2         : int  78524315 78524332 78524325 78524373 78524315 78524315 78524315 78524315 78524541 78524321 ...
#  $ MAPQ_2        : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ CIGAR_2       : chr  "108S82M12D79M" "65M12D86M" "72M12D79M" "24M12D127M" ...
#  $ SEQ_2         : chr  "CACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACGTGCCTTTAAATGAAGATGGTAATTATTACTATTACTGCCACCACCTTTGATAAATTTAGGTCTG"| __truncated__ "CTCCCTCTCCATTAAAATGCAACTTAATTGGTAAAGTTTAGTTGACTAGAGATTTCATTTATTTATTTATTTATTTATTTATTTATTTGTTTGTTTGAGATGAGTCTTGCT"| __truncated__ "CTTTCTTCTCCCTCTCCATTAAAATGCAACTTAATTGGTAAAGTTTAGTTGACTAGAGATTTCATTTATTTATTTATTTATTTATTTATTTATTTGTTTGTTTGAGATGAG"| __truncated__ "TTGACTAGAGATTTCATTTATTTATTTATTTATTTATTTATTTATTTGTTTGTTTGAGATGAGTCTTGCTCTTTTGCCCAAGCTGGAGTGCAGTGGCATGATCTCGGCTCA"| __truncated__ ...
#  $ ALN_N_2       : int  2 2 2 2 1 1 1 1 2 2 ...
#  $ UMI_2         : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ groupIndex    : int  10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 ...
#  $ chrom1        : int  16 16 16 16 16 16 16 16 16 16 ...
#  $ side1         : chr  "R" "R" "R" "R" ...
#  $ pos1          : int  78523150 78523183 78523199 78523164 78523150 78523150 78523150 78523150 78523216 78523183 ...
#  $ chrom2        : int  16 16 16 16 16 16 16 16 16 16 ...
#  $ side2         : chr  "R" "R" "R" "R" ...
#  $ pos2          : int  78524315 78524332 78524325 78524373 78524315 78524315 78524315 78524315 78524541 78524321 ...
#  $ jxnName       : chr  "16:R:78523150,16:R:78524315" "16:R:78523183,16:R:78524332" "16:R:78523199,16:R:78524325" "16:R:78523164,16:R:78524373" ...
#  $ jxnKey        : chr  "1:455840:1" "1:455841:1" "1:455842:1" "2:714401:1" ...
#  $ svIndex       : chr  "10000:1" "10000:1" "10000:1" "10000:1" ...
#  $ sampleSvIndex : chr  "1:10000:1" "1:10000:1" "1:10000:1" "2:10000:1" ...
#  $ AMBIGUOUS     : int  0 1 1 1 0 0 0 0 1 1 ...
#  $ DOWNSAMPLED   : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ N_COLLAPSED   : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ IS_REFERENCE  : num  1 0 0 0 0 0 0 0 0 0 ...
#  $ JXN_SEQ       : chr  "CACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACGTGCCTTTAAATGAAGATGGTAATTATTACTATTACTGCCACCACCTTTGATAAATTTAGGTCTG"| __truncated__ NA NA NA ...
#  $ MERGE_LEN     : int  0 NA NA NA NA NA NA NA NA NA ...
#  $ TARGET_POS_1  : int  2867229 2867293 2867309 2867274 2867229 2867229 2867229 2867229 2867326 2867293 ...
#  $ TARGET_POS_2  : int  2868447 2868454 2868447 2868495 2868407 2868447 2868454 2868509 2868651 2868447 ...
