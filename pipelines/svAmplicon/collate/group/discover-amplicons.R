# iteratively discover amplicon endpoint peaks and associate molecules with them

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("  initializing")
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'DISCOVERY_FILE',
        'ALLOWED_FILE',
        'AMPLICONS_FILE'   
    ),
    integer = c(
        'N_CPU',
        'MATCH_KEEP_DISTANCE',
        'MATCH_REJECT_DISTANCE',
        'MAX_AMPLICONS',
        'MAX_INSERT_SIZE',
        'MIN_MERGE_OVERLAP',
        'READ_LEN'
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
write.table(
    keptNodePairs, 
    file = env$ALLOWED_FILE, 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)
#=====================================================================================

#=====================================================================================
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
write.table(
    amplicons, 
    file = env$AMPLICONS_FILE, 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)
#=====================================================================================
