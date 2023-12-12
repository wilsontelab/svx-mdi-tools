#----------------------------------------------------------------------
# constants relevant to svPore and similar node graph representations of reads
# (as opposed to earlier, single-junction encoding like svCapture)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# node classes
#----------------------------------------------------------------------
svx_nodeClasses <- list(
   'GAP'='0',
   'SPLIT'='1',
   'CLIP'='2'
)
svx_nodeClasses_rev <- names(svx_nodeClasses)

#----------------------------------------------------------------------
# node graph edge types, carried in column "edgeType"
# all _apps_ now adopt this nomenclature for junction/edge codification, for consistency
#----------------------------------------------------------------------
svx_edgeTypes <- list( # edges represent all parts of a read/molecule
    ALIGNMENT     = "A", # the single edge type for a contiguous aligned segment
    TRANSLOCATION = "T", # junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    PROPER        = "P", 
    INTERGENOME   = "G"
)
svx_jxnTypes <- data.table( # junctions are the subset of edges that cross SVs, potentially several per read/molecule
    code = c( # edgeType as used in svPore, others
        "D", # deletion type uses CIGAR-consistent codes
        "I", # insertion type supported in edgeType encoding, maintain CIGAR codes, necessitating adjustment of duplication inversion, etc.
        "U",
        "V",
        "T",
        "?",
        "G"
    ),
    altCode = c( # JXN_TYPE as used in svWGS/svCapture, converted to edgeType==code by SV load functions
        "L", 
        "X", # insertion type not used in this encoding
        "D",
        "I",
        "T",
        "?",
        "G"
    ),
    name = c(
        "Del",
        "Ins",
        "Dup",
        "Inv",
        "Trans",
        "?",
        "IntGen"
    ),
    longName = c(
        "Deletion",
        "Insertion",
        "Duplication",
        "Inversion",
        "Translocation",
        "?",
        "Intergenome"
    ), 
    legend = c(
        TRUE,
        FALSE,
        TRUE,
        TRUE,
        TRUE,
        FALSE,
        TRUE
    ),
    order = 1:7,
    color = c(
        CONSTANTS$plotlyColors$blue,
        CONSTANTS$plotlyColors$black,
        CONSTANTS$plotlyColors$green,
        CONSTANTS$plotlyColors$red,
        CONSTANTS$plotlyColors$orange,
        NA,
        CONSTANTS$plotlyColors$purple
    ),
    lineN = c(
        1,
        -100,
        2,
        3,
        4,
        -100,
        4
    )
)
svx_jxnType_codeToX <- function(codes, col){ # use these helper functions for converting between jxnType columns
    jt <- svx_jxnTypes
    setkey(jt, code)
    jt[codes][[col]]
}
svx_jxnType_altCodeToX <- function(altCodes, col){
    jt <- svx_jxnTypes
    setkey(jt, altCode)
    jt[altCodes][[col]]
}
svx_jxnType_nameToX <- function(names, col){
    jt <- svx_jxnTypes[, name := tolower(name)]
    setkey(jt, name)
    jt[tolower(names)][[col]]
}
svx_jxnType_longNameToX <- function(longNames, col){
    jt <- svx_jxnTypes[, longName := tolower(longName)]
    setkey(jt, longName)
    jt[tolower(longNames)][[col]]
}

#----------------------------------------------------------------------
# junction filters
# all values here should be ~unfiltered, so that apps that don't offer the setting don't apply the filter
#----------------------------------------------------------------------
svx_filterDefaults <- list( 
    Min_SV_Size = 1,
    Max_SV_Size = 0,
    Min_Insert_Size = -10000,
    Max_Insert_Size =  10000,
    Min_Source_Molecules = 1,
    Max_Source_Molecules = 0,
    Min_Sequenced_Molecules = 0,  
    Max_Linked_Junctions = 0, 
    Min_Samples_With_SV = 1,
    Max_Samples_With_SV = 0,
    Unique_To_Sample = "show all SVs",     
    Show_ChrM = "always",  
    Min_Map_Quality = 0,
    Min_Flank_Length = 0,
    SV_Type = c("Del","Dup","Inv"),
    Min_Flank_CNC = 0
)
