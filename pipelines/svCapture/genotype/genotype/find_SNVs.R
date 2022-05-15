# examine all called TT and TA class SVs for the presence of concomitant SNVs/indels

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(parallel)
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOMEX_MODULES_DIR',
        'MODULES_DIR',
        'ACTION_DIR',
        'FIND_PREFIX',
        'GENOTYPE_PREFIX'
    ),
    integer = c(
        'N_CPU'
    ),
    logical = c(
        'DUPLEX_ONLY'
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
sourceScripts(file.path(rUtilDir, 'sequence'), c('general'))
sourceScripts(file.path(env$MODULES_DIR, 'find2', 'find'), c('column_definitions'))
#-------------------------------------------------------------------------------------
matchTypes <- list(
    UNSEQUENCED   = 0, # base was in a read gap or an insertion, nothing to compare
    UNINFORMATIVE = 1, # junction base did not match reference, no haplotype calls existed at this base
    HAPLOTYPE     = 2, # junction base matched at least one source allele (the most common outcome)
    REFERENCE     = 3, # junction base did not match either source allele, but did match the reference
    MISMATCH      = 4  # junction base matched neither a source allele nor the genome reference
)
#=====================================================================================

#=====================================================================================
# load required data
#-------------------------------------------------------------------------------------
message("loading unphased haplotype map")
hapFile <- paste(env$GENOTYPE_PREFIX, 'unphased_haplotypes', 'rds', sep = ".")
hapMap <- readRDS(hapFile)
hapMap[, posKey := paste(CHROM, POS, sep = ":")]
setkey(hapMap, posKey)

message("loading and filtering structural variants")
svsFile <- paste(env$FIND_PREFIX, 'structural_variants', 'rds', sep = ".")
svCalls <- readRDS(svsFile)
svCalls <- svCalls[TARGET_CLASS %in% c("TT", "TA") & JXN_BASES != "*"]

message("loading junction molecules")
molsFile <- paste(env$FIND_PREFIX, 'junction_molecules',  'rds', sep = ".")
jxnMols <- readRDS(molsFile)
if(env$DUPLEX_ONLY) jxnMols <- jxnMols[IS_DUPLEX == 1]
jxnMols[, isOuterClip := NODE_CLASS == nodeClasses$OUTER_CLIP]
unpackNodeNames <- function(nodeNames){ 
    x <- as.data.frame(t(sapply(nodeNames, function(x){
        if(x == "*") c("0", "*", "0") else strsplit(x, '\\W')[[1]]
    })), stringsAsFactors = FALSE)
    x[[1]] <- as.integer(x[[1]])
    x[[3]] <- as.integer(x[[3]])
    x
}
jxnMols[, c('chrom1', 'side1', 'pos1') := unpackNodeNames(NODE_1)]
jxnMols[, c('chrom2', 'side2', 'pos2') := unpackNodeNames(NODE_2)]
#=====================================================================================

#=====================================================================================
# convert an SV and its set of evidence molecules to a base map
# ------------------------------------------------------------------------------------
getJunctionMap <- function(sv, mols){

    # initialize map and dimensions
    refWidth <- nchar(sv$GEN_REF_1)
    faidxPadding <- (refWidth - 1) / 2
    microhomologyLength <- sv[, MICROHOM_LEN]    
    leftRefI  <- faidxPadding + 1
    rightRefI <- leftRefI + 1 - microhomologyLength
    mapWidth  <- refWidth + 1 - microhomologyLength
    nAln <- mols[, .N + sum(!isOuterClip)]
    yOffset <- 1
    baseMap <- matrix(NA, nrow = mapWidth, ncol = nAln)   
    
    # add all alignments to map
    processNode <- function(nodeN, isRef, jxnSide, jxnPos, CIGAR, SEQ, pos){
        if(CIGAR == "*") return(NULL) # the missing node of an outer clip
        aln <- cigarToRefAln(CIGAR, SEQ)           
        nodePosDelta <- if(jxnSide == "L") jxnPos - pos else pos - jxnPos
        mapPos <- if(nodeN == 1){
            leftRefI - nodePosDelta - aln$lengthAlignedOut - aln$leftClip + 1
        } else {
            leftRefI + 1 - microhomologyLength + nodePosDelta - aln$leftClip
        }
        j <- mapPos:(mapPos + aln$lengthOut - 1)
        baseMap[j, yOffset] <<- aln$seq # clip and insertion are lower case, otherwise ACGTN,-,+
        yOffset <<- yOffset + 1 
    }
    mols[, mapply(processNode, 1, IS_REFERENCE, sv$SIDE_1, sv$POS_1, CIGAR_1, SEQ_1, pos1)]
    mols[, mapply(processNode, 2, IS_REFERENCE, sv$SIDE_2, sv$POS_2, CIGAR_2, SEQ_2, pos2)]

    # assemble the final maps
    usedPos <- rle(apply(baseMap[, ], 1, function(v) any(!is.na(v))))$lengths
    usedPos <- (usedPos[1] + 1):(mapWidth - usedPos[length(usedPos)]) # row indices = base positions in SV allele
    firstUsedPos <- usedPos[1]

    # set chrom posiitions, including reverse complement
    isRC1 <- FALSE
    isRC2 <- FALSE
    if(sv$SIDE_1 == sv$SIDE_2){
        if(sv$SIDE_1 == "L") {
            POS_1 <- (sv$POS_1 - faidxPadding):(sv$POS_1 + faidxPadding)
            POS_2 <- (sv$POS_2 + faidxPadding):(sv$POS_2 - faidxPadding)      
            isRC2 <- TRUE
        } else {
            POS_1 <- (sv$POS_1 + faidxPadding):(sv$POS_1 - faidxPadding)
            POS_2 <- (sv$POS_2 - faidxPadding):(sv$POS_2 + faidxPadding)      
            isRC1 <- TRUE
        }
    } else {
        POS_1 <- (sv$POS_1 - faidxPadding):(sv$POS_1 + faidxPadding)
        POS_2 <- (sv$POS_2 - faidxPadding):(sv$POS_2 + faidxPadding)       
    }
    POS_1 <- POS_1[1:mapWidth]
    POS_2 <- rev(rev(POS_2)[1:mapWidth])

    # finish and return the results
    list(
        matrix    = baseMap[usedPos, ], 
        leftRefI  = leftRefI  - firstUsedPos + 1,
        rightRefI = rightRefI - firstUsedPos + 1,
        posKey1   = paste(sv$CHROM_1, POS_1[usedPos], sep = ":"),
        posKey2   = paste(sv$CHROM_2, POS_2[usedPos], sep = ":"),
        isRC1     = isRC1,
        isRC2     = isRC2
    )
}
#=====================================================================================

#=====================================================================================
# convert the base map into a consensus sequence of the junction
# ------------------------------------------------------------------------------------
getJunctionConsensus <- function(matrix){
    apply(matrix, 1, function(x){
        x[x %in% clipBases] <- NA # insertions are also lower case
        x <- x[!is.na(x)]
        if(length(x) == 0) return("N")
        agg <- aggregate(x, list(x), length)
        agg[which.max(agg[[2]]), 1]
    })
}
#=====================================================================================

#=====================================================================================
# compare the junction consensus to the reference and source haplotypes
# ------------------------------------------------------------------------------------
#     JXN     N
# 1:    C  7403
# 2:    T 12717
# 3:    A 12621
# 4:    G  7271
# 5: <NA>   479
# 6:    N    11
# 7:    +     1
# 8:    -    13
compareSequences <- function(REF, HAP1, HAP2, INF, JXN){
    if("N" %in% c(HAP1, HAP2, JXN)) return(matchTypes$UNSEQUENCED)
    if(JXN == HAP1 || # includes both base matches and deleted bases
       JXN == HAP2 ||
       (nchar(HAP1) > 1 && JXN == "+") || # not matching insertion sequences, just presence
       (nchar(HAP2) > 1 && JXN == "+")
    ) return(matchTypes$HAPLOTYPE)
    if(INF == 0) return(matchTypes$UNINFORMATIVE) # thus, can match, but not mismatch, if uninformative
    if(JXN == REF) return(matchTypes$REFERENCE)
    matchTypes$MISMATCH
}
simplifyHaplotype <- Vectorize(function(REF, HAP){
    if(REF == HAP) "." else HAP
})
parseMatches <- Vectorize(function(MATCH){
    switch(
        MATCH + 1,
        "~",
        "~",
        "|",
        ".",
        "X"
    )
})
pasteSequence <- function(bases, isHaplotype){
    bases <- ifelse(is.na(bases),     "~", bases)
    bases <- ifelse(nchar(bases) > 1, "+", bases)
    paste(bases, collapse = "")
}
#=====================================================================================

#=====================================================================================
# search for novel SNVs and indels in each SV junction
#-------------------------------------------------------------------------------------
message("comparing junction sequence to source alleles")
svData <- do.call(rbind, mclapply(svCalls$SV_ID, function(svId){
    sv   <- svCalls[SV_ID == svId]
    mols <- jxnMols[SV_ID == svId]
    if(nrow(mols) == 0) return(NULL) # e.g., if DUPLEX_ONLY is set and is not duplex
    if(mols[, sum(!isOuterClip)] == 0) return(NULL) # or the only duplex molecules were outer clips
    jxnMap <- getJunctionMap(sv, mols)
    consensus <- getJunctionConsensus(jxnMap$matrix)
    pos1 <- 1:jxnMap$leftRefI
    pos2 <- jxnMap$rightRefI:length(jxnMap$posKey2)
    hapMap1 <- hapMap[jxnMap$posKey1[pos1]]
    hapMap2 <- hapMap[jxnMap$posKey2[pos2]]
    hapMap1[, JXN := consensus[pos1]]
    hapMap2[, JXN := consensus[pos2]]
    if(jxnMap$isRC1) hapMap1[, JXN := rc(JXN)]
    if(jxnMap$isRC2) hapMap2[, JXN := rc(JXN)]
    hapMap1[, MATCH := mapply(compareSequences, REF, HAP1, HAP2, INF, JXN)]
    hapMap2[, MATCH := mapply(compareSequences, REF, HAP1, HAP2, INF, JXN)]
    matchType <- max(hapMap1[, max(MATCH)], hapMap1[, max(MATCH)])
    hapMap1[, HAP1 := simplifyHaplotype(REF, HAP1)]
    hapMap1[, HAP2 := simplifyHaplotype(REF, HAP2)]
    hapMap2[, HAP1 := simplifyHaplotype(REF, HAP1)]
    hapMap2[, HAP2 := simplifyHaplotype(REF, HAP2)]
    hapMap1[, MATCH := parseMatches(MATCH)]
    hapMap2[, MATCH := parseMatches(MATCH)]
    data.table(
        SV_ID = svId,
        MATCH_TYPE = matchType,
        INF_1   = pasteSequence(hapMap1$INF),
        REF_1   = pasteSequence(hapMap1$REF),
        HAP1_1  = pasteSequence(hapMap1$HAP1),
        HAP2_1  = pasteSequence(hapMap1$HAP2),
        JXN_1   = pasteSequence(hapMap1$JXN),
        MATCH_1 = pasteSequence(hapMap1$MATCH),
        INF_2   = pasteSequence(hapMap2$INF),
        REF_2   = pasteSequence(hapMap2$REF),
        HAP1_2  = pasteSequence(hapMap2$HAP1),
        HAP2_2  = pasteSequence(hapMap2$HAP2),
        JXN_2   = pasteSequence(hapMap2$JXN),
        MATCH_2 = pasteSequence(hapMap2$MATCH)
    )
}, mc.cores = env$N_CPU)) # 
#=====================================================================================

#=====================================================================================
# print results
#-------------------------------------------------------------------------------------
message("writing SV summary table")
outFile <- paste(env$GENOTYPE_PREFIX, 'haplotype_comparisons', 'gz', sep = ".")
fwrite(
    svData, 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = TRUE, 
    compress = "gzip"
)
outFile <- paste(env$GENOTYPE_PREFIX, 'haplotype_comparisons', 'rds', sep = ".")
saveRDS(
    svData, 
    file = outFile
)
#=====================================================================================
  