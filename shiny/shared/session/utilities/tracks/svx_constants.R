# junction types

svx_edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    PROPER        = "P"
)
svx_jxnTypes <- data.table(
    code = c( # as used in svPore
        "D",
        "I",
        "U",
        "V",
        "T",
        "?"
    ),
    altCode = c( # as used in svWGS/capture
        "L",
        "X", # not used
        "D",
        "I",
        "T",
        "?"
    ),
    name = c(
        "Del",
        "Ins",
        "Dup",
        "Inv",
        "Trans",
        "?"
    ),
    longName = c(
        "Deletion",
        "Insertion",
        "Duplication",
        "Inversion",
        "Translocation",
        "?"
    ), 
    legend = c(
        TRUE,
        FALSE,
        TRUE,
        TRUE,
        TRUE,
        FALSE
    ),
    order = 1:6,
    color = c(
        CONSTANTS$plotlyColors$blue,
        CONSTANTS$plotlyColors$black,
        CONSTANTS$plotlyColors$green,
        CONSTANTS$plotlyColors$red,
        CONSTANTS$plotlyColors$orange,
        NA
    ),
    lineN = c(
        1,
        -100,
        2,
        3,
        4,
        -100
    )
)
svx_jxnType_codeToX <- function(x, col){
    jt <- svx_jxnTypes
    setkey(jt, code)
    jt[x][[col]]
}
svx_jxnType_altCodeToX <- function(x, col){
    jt <- svx_jxnTypes
    setkey(jt, altCode)
    jt[x][[col]]
}
svx_jxnType_nameToX <- function(x, col){
    jt <- svx_jxnTypes
    setkey(jt, name)
    jt[x][[col]]
}

# junction filters
svx_filterDefaults <- list(
    Min_SV_Size = 1,
    Max_SV_Size = 0,
    Min_Insert_Size = -50,
    Max_Insert_Size = 50,
    Min_Samples_With_SV = 1,
    Max_Samples_With_SV = 0,
    Min_Source_Molecules = 1,
    Max_Source_Molecules = 0,
    Min_Sequenced_Molecules = 0,  
    Max_Linked_Junctions = 0, 
    Min_Map_Quality = 0,
    Min_Flank_Length = 0,
    SV_Type = c("Del","Dup","Inv"),  
    Show_ChrM = "never"    
)
svx_nodeClasses <- list(
   'GAP'='0',
   'SPLIT'='1',
   'CLIP'='2'
)
svx_nodeClasses_rev <- names(svx_nodeClasses)
