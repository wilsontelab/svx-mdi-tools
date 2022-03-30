
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
        'FIND_PREFIX',
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
        'MIN_COVERAGE',
        'MIN_MERGE_OVERLAP'
    ),
    double = c(
        'MIN_MERGE_DENSITY'
    )   
))
writeEnvJson(env$FIND_PREFIX) # for use by Shiny app
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
sourceScripts(file.path(env$ACTION_DIR, 'find'), c('column_definitions', 'node_retrieval', 'network', 'analyze_junctions'))
#-------------------------------------------------------------------------------------
# initialize genome
setCanonicalChroms()
loadFaidx(env$SHM_DIR_WRK)
#-------------------------------------------------------------------------------------
# initialize samples
env$SAMPLES   <- strsplit(env$SAMPLES,   "\\s+")[[1]]
SAMPLES <- as.list(1:length(env$SAMPLES))
names(SAMPLES) <- env$SAMPLES
env$MAX_TLENS <- strsplit(env$MAX_TLENS, "\\s+")[[1]]
MAX_MAX_TLEN <- max(env$MAX_TLENS)
MAX_TLENS <- as.list(env$MAX_TLENS)
names(MAX_TLENS) <- env$SAMPLES
faidx_padding <- round(MAX_MAX_TLEN * 1.2, 0) # sufficient to contain any source molecule span
#=====================================================================================

#=====================================================================================
# load and preprocess the called and grouped SVs across all samples
#-------------------------------------------------------------------------------------
jxnMols <- as.data.table(read.table(
    'stdin',
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    ######################
    col.names  = names(find$working),
    colClasses = unname(unlist(find$working))
))
# setkey(jxns, groupIndex)
# groupIndices <- jxns[, unique(groupIndex)]
jxnMols[, ':='(
    AMBIGUOUS = 0, # TRUE (1) if a junction molecule is in more than one SV group
    jxnName   = paste(NODE_1, NODE_2, sep = ","), # the junction edge called by a molecule (could be >1 per molecule)
    jxnKey    = paste(SAMPLES[SAMPLE], MOL_ID, JXN_N, sep = ":") # ID for the source junction edge (each with 2 nodes)
)]
#=====================================================================================

#=====================================================================================
# parallel process continuity groups to find SV matches between junction molecules
# these steps are executed without respect to source sample
#-------------------------------------------------------------------------------------
# compare a gap junction in one molecule to a sequenced junction in another molecule
checkHalfSequenced <- function(delta, side) { 
    if(side == "L") delta >= -env$PURGE_DISTANCE && delta <= env$MAX_MAX_TLEN
               else delta >= -env$MAX_MAX_TLEN   && delta <= env$PURGE_DISTANCE  
}
# determine if edges in different molecules appear to be the same SV on one side of the junction
isJxnMatch <- function(IS_SPLIT_J1, IS_SPLIT_J2, POS_J1, POS_J2, side){
         if( IS_SPLIT_J1 &&  IS_SPLIT_J2) abs(POS_J1 - POS_J2) <= env$PURGE_DISTANCE
    else if(!IS_SPLIT_J1 && !IS_SPLIT_J2) abs(POS_J1 - POS_J2) <= env$MAX_MAX_TLEN
    else if(IS_SPLIT_J1) checkHalfSequenced(POS_J1 - POS_J2, side)
                    else checkHalfSequenced(POS_J2 - POS_J1, side)
}
# determine which jxnMols match a seed junction
addMatches <- function(matches, jxnMols, seedJxn, side1, side2){
    SEED_IS_SPLIT <- jxnMols[seedJxn, NODE_CLASS == nodeClasses$SPLIT]
    IS_SPLIT      <- jxnMols[       , NODE_CLASS == nodeClasses$SPLIT] 
    cbind(matches, 
        mapply(
            isJxnMatch, 
            SEED_IS_SPLIT, 
            IS_SPLIT, 
            jxnMols[seedJxn, POS_1], 
            jxnMols[       , POS_1],
            side1
        ) & mapply(
            isJxnMatch, 
            SEED_IS_SPLIT, 
            IS_SPLIT, 
            jxnMols[seedJxn, POS_2], 
            jxnMols[       , POS_2],
            side2
        )
    )
}
commitSV <- function(jxnMols){
    jxnMols #############################
}

# determine which junction molecules in a continuity group appear to be the same SV
parseContinuityGroup <- function(jxnMols){

    # check if singleton-molecule call or a single-junction set, e.g., multiple identical splits
    if(nrow(jxnMols) == 1 ||
       jxnMols[, length(unique(jxnName)) == 1]) return(commitSV(jxnMols))

    # pick a seed junction, preferring splits > gaps, higher coverage
    setkey(jxnMols, jxnName)
    seedJxn <- jxnMols[order(-NODE_CLASS, -(STRAND_COUNT1 + STRAND_COUNT2)), jxnName[1]]

    # mark matches to seed
    matches <- data.table()
    side1 <- jxnMols[1, strsplit(NODE_1, ":")[[1]][2] ] # sides are the same throughout a continuity group
    side2 <- jxnMols[1, strsplit(NODE_2, ":")[[1]][2] ]
    matches <- addMatches(matches, jxnMols, seedJxn, side1, side2)

    # check if one coherent SV group, some combination of splits and gaps
    if(all(matches)) return(commitSV(jxnMols))

    # pick the next best seed from the unmatched junctions, iterate until all junctions match a seed
    matchedJxns <- apply(matches, 1, any)
    while(!all(matchedJxns)){
        seedJxn <- jxnMols[matchedJxns == FALSE][order(-NODE_CLASS, -(STRAND_COUNT1 + STRAND_COUNT2)), jxnName[1]]
        matches <- addMatches(matches, jxnMols, seedJxn, side1, side2)
        matchedJxns <- apply(matches, 1, any)
    }

    # analyze the SV groups for ambiguous gap-junction molecules consistent with multiple seeds
    # these arise when:
    #   two junctions are closely spaced, such that a gap could have been from either allele source
    #   an extremely large source molecule failed MAX_TLEN matching, such that it became a false seed
    #   a split node has an alignment position error more extreme than PURGE_DISTANCE, another false seed
    nSeedMatches <- apply(matches, 1, sum)
    jxnMols[, AMBIGUOUS := as.integer(nSeedMatches > 1)]

    # commit every seed junction group as an SV; ambiguous molecules will be present multiple times
    do.call(rbind, lapply(1:ncol(matches), function(i) 
        commitSV(jxnMols[matches[, ..i] == TRUE])
    ))

    # # check for a junction match between every pair of SV calls across all samples
    # pairs <- as.data.table(t(combn(jxns[, jxnKey], 2)), stringsAsFactors = FALSE)
    # setnames(pairs, c('jxnKey1', 'jxnKey2')) 
    # pairs <- cbind(
    #     pairs, 
    #     jxns[pairs[,jxnKey1], .(IS_SPLIT_J1 = !is.na(JXN_SEQ), POS_1_J1 = POS_1, POS_2_J1 = POS_2)], 
    #     jxns[pairs[,jxnKey2], .(IS_SPLIT_J2 = !is.na(JXN_SEQ), POS_1_J2 = POS_1, POS_2_J2 = POS_2)]
    # )
    # pairs[, # check each side of the junction independently
    #     collision := mapply(checkCollision, IS_SPLIT_J1, IS_SPLIT_J2, POS_1_J1, POS_1_J2, jxns[1, SIDE_1]) &
    #                  mapply(checkCollision, IS_SPLIT_J1, IS_SPLIT_J2, POS_2_J1, POS_2_J2, jxns[1, SIDE_2])
    # ]
    # if(pairs[, !any(collision)]) return(x)

    # # find and process all of the unique SV groups in the continuity group

    # # if SV matches were found, concatenate the list of other matching SVs for each SV
    # pairs <- rbind(
    #     pairs[collision == TRUE, .(svKey = svKey1, otherSv = svKey2)],
    #     pairs[collision == TRUE, .(svKey = svKey2, otherSv = svKey1)]
    # )
    # pairs <- pairs[, .(
    #     N_MATCHES_   = length(otherSv),
    #     MATCHING_SVS_ = paste(otherSv, collapse=",")
    # ), by = svKey]
    # x[pairs[, svKey], ':='(
    #     N_MATCHES    = pairs[, N_MATCHES_],
    #     MATCHING_SVS = pairs[, MATCHING_SVS_]
    # )]

    # # return the modified SV list
    # x
}
#-------------------------------------------------------------------------------------
# svs <- do.call(rbind, mclapply(groupIndices, markSharedSVs, mc.cores = env$N_CPU))
svs <- jxnMols[,
    parseContinuityGroup(.SD), 
    keyby = groupIndex # analyze each continuity group separately for one or more SVs
]
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
