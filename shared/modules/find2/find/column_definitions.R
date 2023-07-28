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
clipBases <- c("a", "c", "g", "t", "n")

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
compile <- list(nodes = list(
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
compile$junctions <- c(compile$nodes, list(
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
# working file columns (find's extensions of compiled junction molecules)
#----------------------------------------------------------------------
find <- list(working = c(compile$junctions, list(
    'groupIndex'    = 'integer'
)))
find$working2 <- c(find$working, list(
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

#----------------------------------------------------------------------
# finde output file columns
#----------------------------------------------------------------------
find$structural_variants <- list( # FOR REFERENCE ONLY; last columns are sample-specific, not standarized
    'SV_ID'     =  'character',  # this column has a changed order relative to analyze_junctions.R:characterizeSvJunction() due to sample count merge
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
    'FLANK_LEN1'    =  'integer',
    'FLANK_LEN2'    =  'integer',
    'N_CLUSTERED_JUNCTIONS' = 'integer',
    #---------------
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
find$junction_molecules <- c(compile$junctions, list(
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
