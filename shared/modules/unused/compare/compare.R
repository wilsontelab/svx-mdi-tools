
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
        'COMPARE_PREFIX',        
        'GENOME',
        'GENOME_CHROMS'
    ),
    integer = c(
        'N_CPU',
        'SPLIT_ALLOWANCE',
        'MAX_MAX_TLEN'
    )
))
writeEnvJson(env$COMPARE_PREFIX) # for use by compare and/or Shiny app
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn=2) ########################
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(env$MODULES_DIR, 'find', 'find'), c('column_definitions'))
#=====================================================================================

#=====================================================================================
# load and preprocess the called and grouped SVs across all samples
#-------------------------------------------------------------------------------------
svs <- as.data.table(read.table(
    'stdin',
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names  = names(compare$working),
    colClasses = unname(unlist(compare$working))
))
groupIndices <- svs[, unique(groupIndex)]
#=====================================================================================

#=====================================================================================
# parallel process continuity groups to find SV collisions between samples
#-------------------------------------------------------------------------------------
# compare a gap junction in one sample to a sequenced junction in the other sample
checkHalfSequenced <- function(delta, side) { 
    if(side == "L") delta >= -env$SPLIT_ALLOWANCE && delta <= env$MAX_MAX_TLEN
               else delta >= -env$MAX_MAX_TLEN    && delta <= env$SPLIT_ALLOWANCE  
}
# determine if calls in different samples appear to be the same SV on one side of the junction
checkCollision <- function(IS_SEQ_S1, IS_SEQ_S2, POS_S1, POS_S2, side){
         if( IS_SEQ_S1 &&  IS_SEQ_S2) abs(POS_S1 - POS_S2) <= env$SPLIT_ALLOWANCE
    else if(!IS_SEQ_S1 && !IS_SEQ_S2) abs(POS_S1 - POS_S2) <= env$MAX_MAX_TLEN
    else if(IS_SEQ_S1) checkHalfSequenced(POS_S1 - POS_S2, side)
                  else checkHalfSequenced(POS_S2 - POS_S1, side)
}
# determine which SVs in a continuity group appear to be the same SV
markSharedSVs <- function(groupI){

    # if just one SV in a continuity group, commit it as unique to one sample
    x <- svs[groupIndex == groupI] 
    if(nrow(x) == 1) return(x)
    setkey(x, svKey)

    # check for a junction match between every pair of SV calls across all samples
    pairs <- as.data.table(t(combn(x[, svKey], 2)), stringsAsFactors = FALSE)
    setnames(pairs, c('svKey1', 'svKey2')) 
    pairs <- cbind(
        pairs, 
        x[pairs[,svKey1], .(IS_SEQ_S1 = !is.na(JXN_SEQ), POS_1_S1 = POS_1, POS_2_S1 = POS_2)], 
        x[pairs[,svKey2], .(IS_SEQ_S2 = !is.na(JXN_SEQ), POS_1_S2 = POS_1, POS_2_S2 = POS_2)]
    )
    pairs[, # check each side of the junction independently
        collision := mapply(checkCollision, IS_SEQ_S1, IS_SEQ_S2, POS_1_S1, POS_1_S2, x[1, SIDE_1]) &
                     mapply(checkCollision, IS_SEQ_S1, IS_SEQ_S2, POS_2_S1, POS_2_S2, x[1, SIDE_2])
    ]
    if(pairs[, !any(collision)]) return(x)

    # if SV matches were found, concatenate the list of other matching SVs for each SV
    pairs <- rbind(
        pairs[collision == TRUE, .(svKey = svKey1, otherSv = svKey2)],
        pairs[collision == TRUE, .(svKey = svKey2, otherSv = svKey1)]
    )
    pairs <- pairs[, .(
        N_MATCHES_   = length(otherSv),
        MATCHING_SVS_ = paste(otherSv, collapse=",")
    ), by = svKey]
    x[pairs[, svKey], ':='(
        N_MATCHES    = pairs[, N_MATCHES_],
        MATCHING_SVS = pairs[, MATCHING_SVS_]
    )]

    # return the modified SV list
    x
}
#-------------------------------------------------------------------------------------
svs <- do.call(rbind, mclapply(groupIndices, markSharedSVs, mc.cores = env$N_CPU))
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
