# iteratively discover amplicon endpoint peaks and associate molecules with them

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("  initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)    
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'DISCOVERY_FILE',
        'ALLOWED_FILE',
        'AMPLICONS_FILE',
        'GENOME_FASTA',
        'PLOTS_DIR',
        'PLOT_PREFIX'
    ),
    integer = c(
        'N_CPU',
        'MATCH_KEEP_DISTANCE',
        'MATCH_REJECT_DISTANCE',
        'MAX_AMPLICONS',
        'MAX_INSERT_SIZE',
        'MIN_MERGE_OVERLAP',
        'READ_LEN',
        'PRIMER_MATCH_LENGTH'
    ),
    double = c(
        'MIN_FRACTION_ADD_INDEX'
    )   
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
message("  loading node count file")
nodePairs <- fread(env$DISCOVERY_FILE, sep = "\t", header = FALSE)
setnames(nodePairs, c("chrom1","side1","pos1","chrom2","side2","pos2","count"))
nodePairs[, ":="(
    amplicon = 0,
    index = 0
)]

message("  finding index nodes")
firstAmplicon <- nodePairs[1]
nDiscoveredAmplicons <- 0
keptNodePairs <- nodePairs[FALSE]
while(nodePairs[1, count >= firstAmplicon$count * env$MIN_FRACTION_ADD_INDEX] && 
      nDiscoveredAmplicons < env$MAX_AMPLICONS) {
    nDiscoveredAmplicons <- nDiscoveredAmplicons + 1
    nodePairs[1, index := 1]
    idx <- nodePairs[1]
    matchingI <- nodePairs[, # includes the index node itself
        chrom1 == idx$chrom1 & 
        side1  == idx$side1 & 
        pos1 %between% c(idx$pos1 - env$MATCH_REJECT_DISTANCE, 
                         idx$pos1 + env$MATCH_REJECT_DISTANCE) & 
        chrom2 == idx$chrom2 & 
        side2  == idx$side2 & 
        pos2 %between% c(idx$pos2 - env$MATCH_REJECT_DISTANCE, 
                         idx$pos2 + env$MATCH_REJECT_DISTANCE)
    ]
    nodePairs[matchingI, amplicon := nDiscoveredAmplicons]
    keptI <- nodePairs[,
        matchingI & 
        pos1 %between% c(idx$pos1 - env$MATCH_KEEP_DISTANCE, 
                         idx$pos1 + env$MATCH_KEEP_DISTANCE) & 
        pos2 %between% c(idx$pos2 - env$MATCH_KEEP_DISTANCE, 
                         idx$pos2 + env$MATCH_KEEP_DISTANCE) 
    ]
    keptNodePairs <- rbind(keptNodePairs, nodePairs[keptI])
    nodePairs <- nodePairs[-which(matchingI)]
}
#=====================================================================================

#=====================================================================================
message("  analyzing amplicons")
amplicons <- keptNodePairs[index == 1]
amplicons[, proper := ifelse(
    chrom1 != chrom2 | 
    side1 == side2 | 
    side1 == "L" |
    pos2 - pos1 > env$MAX_INSERT_SIZE, 
    "notPossible",
    ifelse(
        pos2 - pos1 < 2 * env$READ_LEN - env$MIN_MERGE_OVERLAP,
        "expectOverlap",
        "expectGaps"
    )
)]
getChromSpan <- function(chrom, pos1, pos2){
    x <- system2("samtools", c("faidx", env$GENOME_FASTA, paste0(chrom, ":", pos1, "-", pos2)), stdout = TRUE)
    x <- paste(x[2:length(x)], collapse = "")
    toupper(x)
}
getPosSide <- function(chrom, side, pos){
    if(side == L) getChromSpan(chrom, pos - env$MAX_INSERT_SIZE + 1, pos)
             else getChromSpan(chrom, pos, pos + env$MAX_INSERT_SIZE -1)
}
amplicons[, ":="(
    ref1 = if(proper == "notPossible") getPosSide(chrom1, side1, pos1) else getChromSpan(chrom1, pos1, pos2),
    ref2 = if(proper == "notPossible") getPosSide(chrom2, side2, pos2) else "*"
), by = amplicon]
setPaddingBases <- function(side, ref, alt){
    if(ref == "*") ref <- alt
    if(env$MATCH_KEEP_DISTANCE < 1) return( list() )
    length <- nchar(ref)
    sapply(1:env$MATCH_KEEP_DISTANCE, function(nBases){
        if(side == "R") substr(ref, 1, nBases)
                   else substr(ref, length - nBases + 1, length)
    })
}
getPrimer <- function(side, ref, alt){
    if(ref == "*") ref <- alt
    length <- nchar(ref)
    if(side == "R") substr(ref, 1, env$PRIMER_MATCH_LENGTH)
               else substr(ref, length - env$PRIMER_MATCH_LENGTH + 1, length)
}
amplicons[, ":="(
    padding1 = list(setPaddingBases(side1, ref1)),
    padding2 = list(setPaddingBases(side2, ref2, ref1)),
    primer1  = getPrimer(side1, ref1),
    primer2  = getPrimer(side2, ref2, ref1) # NOT reverse complemented
), by = amplicon]
write.table(
    amplicons[, .SD, .SDcols = c(
        "amplicon","proper","count",
        "chrom1","side1","pos1","ref1","primer1",
        "chrom2","side2","pos2","ref2","primer2"
    )], 
    file = env$AMPLICONS_FILE, 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)
#=====================================================================================

#=====================================================================================
message("  repairing molecule outer endpoints to match amplicons")
getEndPatch <- Vectorize(function(amplicon_, refPos, side, padding, pos){
    amplicon <- amplicons[amplicon == amplicon_]
    if(amplicon[[refPos]] == pos) return("0")
    delta <- if(amplicon[[side]] == "R") pos - amplicon[[refPos]] else amplicon[[refPos]] - pos
    if(delta < 0) return(as.character(delta))
    amplicon[[padding]][[1]][delta]
})
keptNodePairs[, ":="(
    patch1 = unlist(getEndPatch(amplicon, "pos1", "side1", "padding1", pos1)),
    patch2 = unlist(getEndPatch(amplicon, "pos2", "side2", "padding2", pos2))
)]
write.table(
    keptNodePairs, 
    file = env$ALLOWED_FILE, 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)
#=====================================================================================
