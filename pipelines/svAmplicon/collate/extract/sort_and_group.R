# sort and group molecules based on the SVs they carry

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
        'OUTPUT_DIR',
        'DATA_NAME',
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
# set some utilities and constants
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
nodeClasses <- c(
    GAP   = 0, # SV evidence type codes, i.e. node classes
    SPLIT = 1    
)
#=====================================================================================

#=====================================================================================
message("  loading interim file")
molecules <- fread(env$INTERIM_FILE, sep = "\t", header = FALSE)
columnsIn <- list(
    amplicon    = 1, # name = column order from parse_nodes.pl, value = output column order
    jxns        = 19,
    indels      = 20,
    isRef       = 5,
    nReadPairs  = 6,
    alns        = 18,
    molTypeId   = 2, # was called molId upstream, will be become the type id once an index is selected
    molClass    = 13,
    isMerged    = 3,
    overlap     = 4,
    tLen        = 8,
    nMBases     = 9,
    nDBases     = 10,
    nIBases     = 11,
    avgQual     = 12,
    nAlns       = 14,
    nJxns       = 15,
    jxnsKey     = 16,
    maxIntMapQ  = 17,
    seq1        = 21,
    qual1       = 22,
    seq2        = 23,
    qual2       = 24,
    nMols       = 7
)
columnsOut <- names(columnsIn)[order(unlist(columnsIn))]
setnames(molecules, names(columnsIn)[names(columnsIn) != "nMols"])
#=====================================================================================

#=====================================================================================
message("  sorting and grouping molecule types")
moleculeTypes <- molecules[
    order(amplicon, jxns, indels, # first three columns are for grouping
          -isRef, -nReadPairs)    # last two columns are for selecting the index molecule
][,
    .(
        # amplicon    = amplicon[1], # molecule-level properties
        molTypeId   = molTypeId[1],  # one molecule id selected as the index for the type
        isMerged    = as.logical(getmode(isMerged)),
        overlap     = getmode(overlap),
        isRef       = as.logical(isRef[1]),
        nReadPairs  = sum(nReadPairs), # continue summing all original read pairs that has this type   
        nMols       = .N, # the number of distinct read pair _sequences_ that had this type
        #---------------------------
        tLen        = tLen[1], # base-level properties
        nMBases     = nMBases[1],
        nDBases     = nDBases[1],
        nIBases     = nIBases[1],        
        avgQual     = max(avgQual), # here and below, a molecule type is supported as the best quality that called it
        #---------------------------
        molClass    = molClass[1], # SV characterization
        nAlns       = nAlns[1],
        nJxns       = nJxns[1],
        jxnsKey     = jxnsKey[1],
        maxIntMapQ  = max(maxIntMapQ),
        #---------------------------
        alns        = alns[1],    
        # jxns        = jxns[1],
        # indels      = indels[1],
        #---------------------------
        seq1        = seq1[1],
        qual1       = qual1[1],
        seq2        = seq2[1],
        qual2       = qual2[1]
    ), by = .(amplicon, jxns, indels) # group
][
    , .SD, .SDcols = columnsOut # reorder columns
][
    order(-nReadPairs) # sort by molecule type frequency
]
#=====================================================================================

#=====================================================================================
message("  creating a table of unique SV junctions in one or more molecules")
junctions <- moleculeTypes[, .(
    level2 = c(
        strsplit(jxns, ":::")[[1]], # consider both splits and small CIGAR indels as an SV-supporting junction
        strsplit(indels, ":::")[[1]]
    )
), by = molTypeId][level2 != "*"][, .(
    level1 = strsplit(level2, "::")[[1]]
), by = .(molTypeId, level2)][level1 != "*"][, {
    jxn <- strsplit(level1, ":")[[1]]
    # join(":", $nodeClass, $jxnType, $svSize, ${$innData[READ1]}[_NODE], ${$innData[READ2]}[_NODE], $overlap, $jxnBases);
    .(
        jxnKey = paste(jxn[2:5], collapse = ":"),
        overlap = as.integer(if(jxn[6] == "NA") NA else jxn[6]),
        jxnBases= as.character(if(jxn[7] == "NA") NA else jxn[7]),
        callType = if(jxn[6] == "NA") "CIGAR" else names(nodeClasses)[nodeClasses == as.integer(jxn[1])]
    )
}, by = .(molTypeId, level1)][jxnKey != "*", .SD, .SDcols = c("molTypeId", "jxnKey","overlap","jxnBases","callType")] # a long-form list correlating all molecules to all junctions

# merge some molecule back fields back into the junction table 
junctions <- merge(
    junctions,
    moleculeTypes[, .SD, .SDcols = c("molTypeId","amplicon","nReadPairs","nMols","isMerged","tLen")],
    by = "molTypeId",
    all.x = TRUE
)

# collapse to lists in both directions, i.e., collapse the jxnKeys within each molTypeId and the molTypeIds that match each jxnKey
moleculeTypes <- merge(
    moleculeTypes, 
    junctions[, .(jxnKeys = list(jxnKey)), by = "molTypeId"], 
    by = "molTypeId",
    all.x = TRUE
)
junctions <- junctions[, 
    .(
        amplicon = amplicon[1],        
        molTypeIds = list(molTypeId),        
        nMolTypes = .N, 
        nMols = sum(nMols),
        nReadPairs = sum(nReadPairs),
        hasMerged =any(isMerged),
        tLens = paste(sort(unique(tLen)), collapse = ","),
        overlap  = if(any(!is.na(overlap)))  overlap[is.na(overlap)][1]   else as.integer(NA),
        jxnBases = if(any(!is.na(jxnBases))) jxnBases[is.na(jxnBases)][1] else as.character(NA),
        callTypes = paste(sort(unique(callType)), collapse = ",")
    ), 
    by = "jxnKey"
]

# split jxnKeys into regular fields for final table displays
junctions <- cbind(
    junctions, 
    t(as.data.table(strsplit(junctions$jxnKey, ":")))
)
junctions <- cbind(
    junctions[, 1:13], 
    t(as.data.table(strsplit(unlist(junctions[, 14]), "/"))),
    t(as.data.table(strsplit(unlist(junctions[, 15]), "/")))
)
setnames(junctions, c(
    "jxnKey",
    "amplicon","molTypeIds","nMolTypes","nMols","nReadPairs","hasMerged","tLens","overlap","jxnBases","callTypes",
    "svType","svSize",
    "chrom1","side1","pos1","chrom2","side2","pos2"
))
junctions[, ":="(
    jxnId  = 1:.N,
    svSize = as.integer(svSize),
    pos1   = as.integer(pos1),
    pos2   = as.integer(pos2)
)]
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
