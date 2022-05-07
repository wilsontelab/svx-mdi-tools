#----------------------------------------------------------------------
# ../column_definitions.R defines the expected extract and find column formats
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# types
#----------------------------------------------------------------------
SVX <- list()
SVX$nodeClasses <- list(
    GAP         = 0, # SV evidence type codes, i.e., node classes
    SPLIT       = 1,
    OUTER_CLIP  = 2
)
SVX$nodeClassColors <- list(
    Gap   = CONSTANTS$plotlyColors$orange,
    Split = CONSTANTS$plotlyColors$blue,
    Clip  = CONSTANTS$plotlyColors$green
)
SVX$getMolColors <- function(x){
    classes <- SVX$nodeClasses
    cols    <- SVX$nodeClassColors
    ifelse(
        x == classes$SPLIT, 
        cols$Split, 
        ifelse(x == classes$GAP, 
            cols$Gap, 
            cols$Clip
        )
    )
}
SVX$jxnTypes <- data.table(
    code = c(
        "D",
        "L",
        "I",
        "T",
        "P",
        "?"
    ),
    name = c(
        "Dup",
        "Del",
        "Inv",
        "Trans",
        "Prop",
        "?"
    ),
    color = c(
        CONSTANTS$plotlyColors$green,
        CONSTANTS$plotlyColors$blue,
        CONSTANTS$plotlyColors$red,
        CONSTANTS$plotlyColors$orange,
        NA,
        NA
    )
)
SVX$targetClasses <- list(
    "both on target" = c("TT", "TA", "AA"),
    "one on target"  = c("TT", "TA", "AA", "tt", "ta", "aa"),
    "all SVs"        = c("TT", "TA", "AA", "tt", "ta", "aa", "t-", "a-", "--")
)

#----------------------------------------------------------------------
# compile file columns
#----------------------------------------------------------------------
unpackNodeNames <- function(nodeNames){ 
    x <- as.data.frame(t(sapply(nodeNames, function(x){
        if(x == "*") c("0", "*", "0") else strsplit(x, '\\W')[[1]]
    })), stringsAsFactors = FALSE)
    x[[1]] <- as.integer(x[[1]])
    x[[3]] <- as.integer(x[[3]])
    x
}
SVX$compile <- list(nodes = list(
    'NODE_1'        = 'character', # node-level data
    'CLIP_LEN_1'    = 'integer',
    'CLIP_SEQ_1'    = 'character',
    #---------------
    'FLAG_1'        = 'integer', # alignment-level data
    'POS_1'         = 'integer',
    'MAPQ_1'        = 'integer',
    'CIGAR_1'       = 'character',
    'SEQ_1'         = 'character',
    'ALN_N_1'       = 'integer',        
    #---------------
    'UMI_1'         = 'integer',
    #===============        
    'NODE_CLASS'    = 'integer',
    #---------------
    'JXN_TYPE'      = 'character', # edge/junction-level data
    'JXN_N'         = 'integer',
    #---------------
    'MOL_ID'        = 'integer', # molecule-level data 
    'IS_MERGED'     = 'integer',
    'IS_DUPLEX'     = 'integer',
    'STRAND_COUNT1' = 'integer',
    'STRAND_COUNT2' = 'integer',
    'MOL_CLASS'     = 'character',
    'MOL_STRAND'    = 'integer',
    'IS_OUTER_CLIP1'= 'integer',
    'IS_OUTER_CLIP2'= 'integer',
    'TARGET_CLASS'  = 'character',
    'SHARED_PROPER' = 'integer',
    #---------------
    'OUT_POS1'      = 'integer',
    'OUT_POS2'      = 'integer'   
))
SVX$compile$junctions <- c(SVX$compile$nodes, list(
    'SAMPLE'        = 'character',
    #===============   
    'NODE_2'        = 'character', # node-level data
    'CLIP_LEN_2'    = 'integer',
    'CLIP_SEQ_2'    = 'character',
    #---------------
    'FLAG_2'        = 'integer', # alignment-level data
    'POS_2'         = 'integer',
    'MAPQ_2'        = 'integer',
    'CIGAR_2'       = 'character',
    'SEQ_2'         = 'character',
    'ALN_N_2'       = 'integer',        
    #---------------
    'UMI_2'         = 'integer'
))

#----------------------------------------------------------------------
# find output file columns
#----------------------------------------------------------------------
SVX$find <- list()
SVX$find$structural_variants <- list( # FOR REFERENCE ONLY; last columns are sample-specific, not standarized
    'SV_ID'     =  'character',  # this column has a changed order relative to analyze_junctions.R:characterizeSvJunction() due to sample count merge # nolint
    #---------------
    'MAPQ_1'    =  'character', # order of these columns generally mirrors junction_molecules
    'UMI_1'     =  'character',  # see analyze_junctions.R:characterizeSvJunction() for the exact content of each column
    #---------------
    'N_TOTAL'   =  'integer',
    'N_GAPS'    =  'integer',
    'N_SPLITS'  =  'integer',
    'N_OUTER_CLIPS' =  'integer',
    #---------------
    'JXN_TYPE'  =  'character',
    #---------------
    'N_DUPLEX'  =  'integer',
    'N_DUPLEX_GS'   =  'integer',
    'STRAND_COUNT'  =  'integer',
    'STRAND_COUNT_GS'   =  'integer',
    'STRAND_COUNT1' =  'integer',
    'STRAND_COUNT2' =  'integer',
    'TARGET_CLASS'  =  'character',
    'SHARED_PROPER' =  'double',
    'SHARED_PROPER_GS'  =  'double',
    #---------------
    'SAMPLES'   =  'character',
    'N_SAMPLES' =  'integer',
    #---------------
    'MAPQ_2'    =  'character',
    'UMI_2'     =  'character',
    #---------------
    'CHROM_1'   =  'character',
    'SIDE_1'    =  'character',
    'POS_1'     =  'integer',
    'CHROM_2'   =  'character',
    'SIDE_2'    =  'character',
    'POS_2'     =  'integer',
    #---------------
    'JUNCTION_NAME'     =  'character',
    'JUNCTION_NAMES' =  'character',
    #---------------
    'N_AMBIGUOUS'   =  'integer',
    'N_DOWNSAMPLED' =  'integer',
    'N_COLLAPSED'   =  'integer',
    #---------------
    'JXN_SEQ'   =  'character',
    'MERGE_LEN' =  'integer',
    #---------------
    'MICROHOM_LEN'  =  'integer',
    'JXN_BASES'     =  'character',
    'SV_SIZE'       =  'integer',
    #---------------
    'GEN_REF_1' =  'character',
    'GEN_REF_2' =  'character',
    'GEN_COV_1' =  'character',
    'GEN_COV_2' =  'character',
    #---------------
    'TARGET_REGION' =  'character',
    'TARGET_POS_1'  =  'integer',
    'TARGET_POS_2'  =  'integer'
    #---------------
    # plus one additional integer column per sample with N_TOTAL for that sample for each SV 
)
SVX$find$junction_molecules <- c(SVX$compile$junctions, list(
    'SV_ID'         =  'character', # for relating this table to find$structural_variants
    #---------------
    'AMBIGUOUS'     = 'integer', # added to compile$junctions by find, per molecule
    'DOWNSAMPLED'   = 'integer', 
    'N_COLLAPSED'   = 'integer',
    'IS_REFERENCE'  = 'integer',
    #---------------
    'TARGET_POS_1'  = 'integer',
    'TARGET_POS_2'  = 'integer'
))
