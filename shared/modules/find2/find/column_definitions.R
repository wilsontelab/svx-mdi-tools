#----------------------------------------------------------------------
# ../column_definitions.R defines the expected extract and find column formats
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# types
#----------------------------------------------------------------------
nodeClasses <- list(
    GAP           = 0, # SV evidence type codes, i.e., node classes
    SPLIT         = 1,
    OUTER_CLIP    = 2,
    RECONSTRUCTED = 3
)
junctionClasses <- list( # i.e., when two nodes have been aggregated
    GAP   = paste(nodeClasses$GAP,   nodeClasses$GAP,   sep = ","),
    SPLIT = paste(nodeClasses$SPLIT, nodeClasses$SPLIT, sep = ","),
    OUTER_CLIP = as.character(nodeClasses$OUTER_CLIP)
)
junctionTypes <- list(
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "I",
    DUPLICATION   = "D",
    DELETION      = "L",
    UNKNOWN       = "?",
    PROPER        = 'P'    
)

#----------------------------------------------------------------------
# compile file columns
#----------------------------------------------------------------------
unpackNodeNames <- function(nodeNames){ 
    x <- as.data.frame(
        matrix(unlist(strsplit(nodeNames, '\\W')), ncol = 3, byrow = TRUE),
        stringsAsFactors = FALSE
    )
    x[[1]] <- as.integer(x[[1]])
    x[[3]] <- as.integer(x[[3]])
    x
}
compile <- list(
    nodes = list(
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
    )
)
compile$junctions <- c(compile$nodes, 
    list(
        #---------------
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
    )  
)
compile$working <- c(compile$junctions, list(
    'groupIndex'    = 'integer'
))
compile$working2 <- c(compile$working, list(
    'chrom1'        = 'integer',
    'side1'         = 'character',
    'pos1'          = 'integer',
    'chrom2'        = 'integer',
    'side2'         = 'character',
    'pos2'          = 'integer',    
    'jxnName'       = 'character',
    'jxnKey'        = 'character', 
    'svIndex'       = 'character',
    'sampleSvIndex' = 'character',    
    'AMBIGUOUS'     = 'integer',
    'DOWNSAMPLED'   = 'integer',
    'N_COLLAPSED'   = 'integer',
    'IS_REFERENCE'  = 'integer'
))

# #----------------------------------------------------------------------
# # find file columns
# #----------------------------------------------------------------------
# find <- list(
#     structural_variants = list( 
#         SV_ID           = "integer",   # SV identifier
#         JUNCTION_NAME   = "character", # the query junction
#         MATCHING_NAMES  = "character", # e.g. gaps that match a split query (should NOT be present as a JUNCTION_NAME in any row) # nolint
#         #-------------
#         TARGET_CLASS    = "character",
#         JXN_TYPE        = "character",
#         #-------------
#         CHROM_1         = "character", # node-level data
#         SIDE_1          = "character",
#         POS_1           = "integer",
#         CHROM_2         = "character",
#         SIDE_2          = "character",
#         POS_2           = "integer",
#         #-------------
#         JXN_SEQ         = "character", # joint (re)construction data
#         MERGE_LEN       = "integer",
#         FAIDX_PADDING   = "integer",
#         GEN_REF_1       = "character",
#         GEN_REF_2       = "character",
#         #-------------
#         MICROHOM_LEN    = "integer",   # joint data
#         MICROHOM_MATCH  = "integer",
#         JXN_BASES       = "character",
#         SV_SIZE         = "integer",     
#         #-------------
#         N_TOTAL = "integer", # evidence counts
#         N_SPLITS = "integer",
#         N_GAPS = "integer",
#         N_OUTER_CLIPS = "integer",
#         N_DUPLEX = "integer", # splits + gaps only
#         N_DUPLEX_ALL = "integer", # including clips
#         NET_STRAND_COUNT = "integer", # splits + gaps only
#         NET_STRAND_COUNT_ALL = "integer", # including clips
#         N_SHARED_PROPER = "integer", # splits + gaps only
#         N_SHARED_PROPER_ALL = "integer", # including clips
#         #------------- 
#         IS_MERGED = "character", # lists across all supporting molecules
#         SEQ_LEN = "character", # is TLEN iff IS_MERGED == TRUE, otherwise is first READ_LEN (e.g. 151)
#         UMI = 'character',
#         MAPQ = 'character',
#         #-------------
#         CHUNK_OFFSET    = "integer",   # index into the all_nodes evidence file
#         CHUNK_SIZE      = "integer",
#         #-------------
#         SAMPLE          = "character" # include sample name for subsequent concatenation during inter-sample comparison
#     ),


# # Classes ‘data.table’ and 'data.frame':  34239 obs. of  48 variables:
# #  $ SV_ID                     : chr  "10000:1" "10000:2" "10000:3" "10000:4" ...
# #  $ MAPQ_1                    : chr  "60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60" "27" "60,60" "60,60" ...
# #  $ UMI_1                     : chr  "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1" "1" "1,1" "1,1" ...
# #  $ N_TOTAL                   : int  20 1 2 2 11 1 19 1 1 1 ...
# #  $ N_GAPS                    : int  10 0 0 1 11 1 19 1 0 0 ...
# #  $ N_SPLITS                  : int  10 1 2 1 0 0 0 0 1 1 ...
# #  $ N_OUTER_CLIPS             : int  0 0 0 0 0 0 0 0 0 0 ...
# #  $ JXN_TYPE                  : chr  "I" "I" "I" "I" ...
# #  $ N_DUPLEX                  : int  0 0 0 0 0 0 0 0 0 0 ...
# #  $ N_DUPLEX_GS               : int  0 0 0 0 0 0 0 0 0 0 ...
# #  $ STRAND_COUNT              : int  157 6 3 11 64 9 118 1 12 1 ...
# #  $ STRAND_COUNT_GS           : int  157 6 3 11 64 9 118 1 12 1 ...
# #  $ STRAND_COUNT1             : int  104 0 1 9 3 0 68 1 0 0 ...
# #  $ STRAND_COUNT2             : int  53 6 2 2 61 9 50 0 12 1 ...
# #  $ TARGET_CLASS              : chr  "TT" "TT" "TT" "TT" ...
# #  $ SHARED_PROPER             : num  1.95 2 2 2 1.82 ...
# #  $ SHARED_PROPER_GS          : num  1.95 2 2 2 1.82 ...
# #  $ SAMPLES                   : chr  "RO3306_0.2APH_Colch_2APH_M,RO3306_0.2APH_Colch_M,RO3306_0.2APH_G2,RO3306_Colch_2APH_M,RO3306_Colch_M,RO3306_G2" "RO3306_0.2APH_G2" "RO3306_Colch_M,RO3306_G2" "RO3306_Colch_M" ...
# #  $ N_SAMPLES                 : int  6 1 2 1 6 1 6 1 1 1 ...
# #  $ MAPQ_2                    : chr  "60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60" "60" "60,45" "60,56" ...
# #  $ UMI_2                     : chr  "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1" "1" "1,1" "1,1" ...
# #  $ CHROM_1                   : chr  "chr16" "chr16" "chr16" "chr16" ...
# #  $ SIDE_1                    : chr  "R" "R" "R" "R" ...
# #  $ POS_1                     : int  78523150 78522944 78522227 78515238 78522947 78519816 78522947 78520486 78514593 78495146 ...
# #  $ CHROM_2                   : chr  "chr16" "chr16" "chr16" "chr16" ...
# #  $ SIDE_2                    : chr  "R" "R" "R" "R" ...
# #  $ POS_2                     : int  78524315 78525450 78525145 78524423 78524737 78524307 78524737 78524757 78526788 78540241 ...
# #  $ JUNCTION_NAME             : chr  "16:R:78523150,16:R:78524315" "16:R:78522944,16:R:78525450" "16:R:78522227,16:R:78525145" "16:R:78515238,16:R:78524423" ...
# #  $ JUNCTION_NAMES            : chr  "16:R:78523150,16:R:78524315::16:R:78523183,16:R:78524332::16:R:78523199,16:R:78524325::16:R:78523164,16:R:78524"| __truncated__ "16:R:78522944,16:R:78525450" "16:R:78522227,16:R:78525145" "16:R:78515467,16:R:78525010::16:R:78515238,16:R:78524423" ...
# #  $ N_AMBIGUOUS               : int  10 0 0 0 7 0 16 0 0 0 ...
# #  $ N_DOWNSAMPLED             : int  0 0 0 0 0 0 0 0 0 0 ...
# #  $ N_COLLAPSED               : int  0 0 0 0 1 0 0 0 0 0 ...
# #  $ JXN_SEQ                   : chr  "CACTATGAAAAATTAATATTTTATTTAGGTCCTTTGGTTTAGATACGTGCCTTTAAATGAAGATGGTAATTATTACTATTACTGCCACCACCTTTGATAAATTTAGGTCTG"| __truncated__ "GCAACCTCCGCCTCCCAGCTTCAAACAATTCTCCTGCCTCAGCCTCCCAGGTAGCTGGGATTACAGGCATTAGNCACCACACCCGGCCCAAATTTTCTTATAAAGTCAAAT"| __truncated__ "TTCAATTAAGCAAAATGGTGTCTTGATTTCACTAGTTTTTTTTTTAATCCTAACCCTTATTTTTCAACTCCTTCTAGGTTTCATGAGCCACTATGCCCAGCCAAATTTTGT"| __truncated__ "GTCTCATCTCAGTGTGCCTCTCAAACACAGATATGCCATTTAATTCCTTAAGTGTTAGTTAAGTAACTACGATGAGCTTAGCTGTACGCAGAGTCTGAAATTACATCTTCA"| __truncated__ ...
# #  $ MERGE_LEN                 : int  0 0 0 0 0 0 0 0 0 0 ...
# #  $ MICROHOM_LEN              : int  13 17 7 18 0 0 0 0 12 9 ...
# #  $ JXN_BASES                 : chr  "" "" "" "" ...
# #  $ SV_SIZE                   : num  1165 2506 2918 9185 1790 ...
# #  $ GEN_REF_1                 : chr  "TGATGGAAGGCTTTAAAAAAAAAAAAAAAAAACTTCAGAAGCTGATTTGAGATTTATTCATCTTTGGTGTACTTTCTTTACTAAGGCAGGTTGTAGACTCATGAAACCTAG"| __truncated__ "TGGGTACTGGAGATGGGGACAGTGAGTATTTGAAGCTCACAGCAGGCCCATTTTGCAAGTTTTCAGGAAACTGAAATTCAGAGAGATTAAGTATTTTGCCATGAATCACTT"| __truncated__ "GATATTATTTGTAGTAATTTATTTAAAGTAATTAGAGACGGGATCTTGCTGTGTTGCCCAGGCTGGTCTCTAACTCCTGGCCTCAAGTGATCCTCCCACCTCAGCCTCCCA"| __truncated__ "AGCAGATTTTAATATCTTTCCTCCCACACTTGCTTTGACTCCCTCCTTTAGTAATTAGATAATTATGTATAATAAGCACATTAATTATGTGCAACAGCTGCCATTGAGTGT"| __truncated__ ...
# #  $ GEN_REF_2                 : chr  "TTCAAGGTAAGGTAGACAGCTGTAAACTGGATCCCCTGCTGTACCCAGCAAATGCTGAGAGCCCTGCTTCCCCTCTTGTATTCAAATGGCTCTATAAGAGAGGAAGTAGTA"| __truncated__ "TTGTTTGTTTGAGATGAGTCTTGCTCTTTTGCCCAAGCTGGAGTGCAGTGGCATGATCTCGGCTCACTGCAGCCTCTGCCTCCCGTGTCCAAGTGATTCTTCTACCTCAGC"| __truncated__ "GTAATTCTCTGATGTGTTTGGGTTGATTCCCTTCTTTCTGCTCTTTTTCTCTTTATAAAGTAAACATTTTTATTGCAAAGGTACTGCTATTGTTAGTTTCTAAAGGCACTT"| __truncated__ "GTACAGCTATCAGTAACAATGCCATCCGGTAATTACCGAATCATAAATGTTTGGGGAGCTGAAGGTCACGGGCTTAGTAGCAATTCTAGTCTGATGTTTGATTCCTTGTGA"| __truncated__ ...
# #  $ TARGET_REGION             : chr  "WWOX" "WWOX" "WWOX" "WWOX" ...
# #  $ TARGET_POS_1              : int  2867110 2866904 2866187 2859198 2866907 2863776 2866907 2864446 2858553 2839106 ...
# #  $ TARGET_POS_2              : int  2868275 2869410 2869105 2868383 2868697 2868267 2868697 2868717 2870748 2884201 ...
# #  $ RO3306_0.2APH_Colch_2APH_M: int  3 0 0 0 1 1 3 0 0 0 ...
# #  $ RO3306_0.2APH_Colch_M     : int  6 0 0 0 2 0 3 0 0 0 ...
# #  $ RO3306_0.2APH_G2          : int  2 1 0 0 1 0 2 0 0 1 ...
# #  $ RO3306_Colch_2APH_M       : int  2 0 0 0 3 0 5 0 0 0 ...
# #  $ RO3306_Colch_M            : int  2 0 1 2 1 0 1 0 1 0 ...
# #  $ RO3306_G2                 : int  5 0 1 0 3 0 5 1 0 0 ...

# )
