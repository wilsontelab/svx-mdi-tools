#----------------------------------------------------------------------
# ../column_definitions.R defines the expected extract and find column formats
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# types
#----------------------------------------------------------------------
edgeTypes <- list(
    junction    = 'J',
    alignment   = 'A',
    termination = 'T'
)
nodeClasses <- list(
    GAP           = 0, # SV evidence type codes, i.e., node classes
    SPLIT         = 1,
    OUTER_CLIP    = 2,
    RECONSTRUCTED = 3
)
edgeClasses <- list(
    sequenced = nodeClasses$SPLIT,
    gap       = nodeClasses$GAP
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
compile <- list(
    edges = list(
        'NODE_CLASS'     = 'integer',
        'NODE1'          = 'character',
        'NODE2'          = 'character',
        'JXN_TYPE'       = 'character',
        'TARGET_CLASS'   = 'character',
        'MOL_IDS'        = 'character',
        'COUNT'          = 'integer',
        'COUNT_DISTINCT' = 'integer' 
    ),
    nodes = list(
        'NODE'           = 'character', # node-level data
        'CLIP_LEN'       = 'integer',
        'CLIP_SEQ'       = 'character',
        'NODE_CLASS'     = 'integer',
        'PARTNER'        = 'character',
        #---------------
        'JXN_TYPE'       = 'character', # edge/junction-level data
        'JXN_N'          = 'integer',
        #---------------
        'FLAG'           = 'integer', # alignment-level data
        'POS'            = 'integer',
        'MAPQ'           = 'integer',
        'CIGAR'          = 'character',
        'SEQ'            = 'character',
        'ALN_N'          = 'integer',
        #---------------
        'MOL_ID'         = 'integer', # molecule-level data 
        'UMI'            = 'integer',
        'IS_MERGED'      = 'integer',
        'IS_DUPLEX'      = 'integer',
        'STRAND_COUNT1'  = 'integer',
        'STRAND_COUNT2'  = 'integer',
        'MOL_CLASS'      = 'character',
        'MOL_STRAND'     = 'integer',
        'IS_OUTER_CLIP1' = 'integer',
        'IS_OUTER_CLIP2' = 'integer',
        'TARGET_CLASS'   = 'character',
        'SHARED_PROPER'  = 'integer',
        'OUT_POS1'       = 'integer',
        'OUT_POS2'       = 'integer'
    )
)

#----------------------------------------------------------------------
# find file columns
#----------------------------------------------------------------------
find <- list(
    structural_variants = list( 
        SV_ID           = "integer",   # SV identifier
        JUNCTION_NAME   = "character", # the query junction
        MATCHING_NAMES  = "character", # e.g. gaps that match a split query (should NOT be present as a JUNCTION_NAME in any row) # nolint
        #-------------
        TARGET_CLASS    = "character",
        JXN_TYPE        = "character",
        #-------------
        CHROM_1         = "character", # node-level data
        SIDE_1          = "character",
        POS_1           = "integer",
        CHROM_2         = "character",
        SIDE_2          = "character",
        POS_2           = "integer",
        #-------------
        JXN_SEQ         = "character", # joint (re)construction data
        MERGE_LEN       = "integer",
        FAIDX_PADDING   = "integer",
        GEN_REF_1       = "character",
        GEN_REF_2       = "character",
        #-------------
        MICROHOM_LEN    = "integer",   # joint data
        MICROHOM_MATCH  = "integer",
        JXN_BASES       = "character",
        SV_SIZE         = "integer",     
        #-------------
        N_TOTAL = "integer", # evidence counts
        N_SPLITS = "integer",
        N_GAPS = "integer",
        N_OUTER_CLIPS = "integer",
        N_DUPLEX = "integer", # splits + gaps only
        N_DUPLEX_ALL = "integer", # including clips
        NET_STRAND_COUNT = "integer", # splits + gaps only
        NET_STRAND_COUNT_ALL = "integer", # including clips
        N_SHARED_PROPER = "integer", # splits + gaps only
        N_SHARED_PROPER_ALL = "integer", # including clips
        #------------- 
        IS_MERGED = "character", # lists across all supporting molecules
        SEQ_LEN = "character", # is TLEN iff IS_MERGED == TRUE, otherwise is first READ_LEN (e.g. 151)
        UMI = 'character',
        MAPQ = 'character',
        #-------------
        CHUNK_OFFSET    = "integer",   # index into the all_nodes evidence file
        CHUNK_SIZE      = "integer",
        #-------------
        SAMPLE          = "character" # include sample name for subsequent concatenation during inter-sample comparison
    ),
    all_nodes = list(
        # logicals here are boolean integers, so they print as 0/1 (for a bit of file space compaction)
        SV_ID            = "integer", # numeric identifier for the SV junction
        NODE_N           = "integer", # which SV junction node this matched
        IS_SEED_NODE     = "integer", # node is an exact match to the the seed/query junction
        IS_REF_NODE      = "integer", # this node was one of two used to characterize the junction
        IS_RC            = "integer", # this SEQ was rc'ed as it was an inversion outer clip
        IS_REPEAT        = "integer", # set in a later step to TRUE/1 if node was already claimed by a lower numbered SV
        N_COLLAPSED      = "integer"  # number of other initial molecules that were merged into this one as presumed duplicates # nolint
    )
)
find$all_nodes <- c(find$all_nodes, compile$nodes)

# the node columns that uniquely identify a specific node-molecule instance
molNodeCols <- c('NODE', 'ALN_N', 'MOL_ID')

#----------------------------------------------------------------------
# compare file columns
#----------------------------------------------------------------------
compare <- list(
    additional = list(
        N_MATCHES    = "integer",  # total number of other SVs matching this SV
        MATCHING_SVS = "character" # comma-delimited list of other SAMPLE:SV_ID that match this SV; '.' if none
    ),
    working = list( # used by compare.R
        svKey      = "character", # SAMPLE:SV_ID
        groupIndex = "integer"    # continuity group marked by makeGroups.pl
    )
)
compare$structural_variants <- c(find$structural_variants, compare$additional)
compare$working <- c(find$structural_variants, compare$additional, compare$working)
