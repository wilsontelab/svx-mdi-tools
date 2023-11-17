# junction types

svx_edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    PROPER        = "P", 
    INTERGENOME   = "G"
)
svx_jxnTypes <- data.table(
    code = c( # as used in svPore
        "D",
        "I",
        "U",
        "V",
        "T",
        "?",
        "G"
    ),
    altCode = c( # as used in svWGS/capture
        "L",
        "X", # not used
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
svx_filterDefaults <- list( # all values here should be ~unfiltered, so that apps that don't offer the setting don't apply the filter
    Min_SV_Size = 1,
    Max_SV_Size = 0,
    Min_Insert_Size = -50,
    Max_Insert_Size = 50,
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
svx_nodeClasses <- list(
   'GAP'='0',
   'SPLIT'='1',
   'CLIP'='2'
)
svx_nodeClasses_rev <- names(svx_nodeClasses)
