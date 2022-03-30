
#----------------------------------------------------------------------
# ../column_definitions.R defines the expected extract and find column formats
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# standard SAM columns
#----------------------------------------------------------------------
SAM_columns <- list(
    'QNAME' = 'character',            # molecule-level information
    'FLAG' = 'integer',
    'RNAME' = 'character',  
    'POS' = 'integer',     
    'MAPQ' = 'integer',
    'CIGAR' = 'character',
    'RNEXT' = 'character',
    'PNEXT' = 'integer',
    'TLEN' = 'integer',
    'SEQ' = 'character',
    'QUAL' = 'character'
)

#----------------------------------------------------------------------
# extract, i.e. molecule, columns
#----------------------------------------------------------------------
extract <- list(
    insertSizeDistribution = list(
        SIZE     = 'integer',
        FREQ     = 'numeric',
        CUM_FREQ = 'numeric'
    ),
    junctions = list(
        'MOL_ID' = 'integer',            # molecule-level information
        'UMI1' = 'integer',
        'UMI2' = 'integer',  
        'IS_DUPLEX' = 'integer',     
        'STRAND_COUNT1' = 'integer',
        'STRAND_COUNT2' = 'integer',
        'IS_MERGED' = 'integer',
        'TARGET_CLASS' = 'character', # target class of the molecule outer positions, not necessrily the SV
        'MOL_CLASS' = 'character',
        'SHARED_PROPER' = 'integer',
        #-------------
        'JXN_CLASS' = 'character',         # junction-level information (class = GAP/SPLIT, type=TDIL)
        'JXN_TYPE' = 'character',
        'WAS_CANONICAL' = 'integer',
        #-------------
        'FLAG1' = 'integer',            # alignment level information
        'RNAME1' = 'character',           # from SAM
        'POS1' = 'integer',
        'MAPQ1' = 'integer',
        'SEQ1' = 'character',
        #-------------
        'OUT_CLIP1' = 'integer',        # outer endpoint data
        'OUT_POS1' = 'integer',
        'OUT_SIDE1' = 'character',
        'OUT_SEQ1' = 'character',
        #-------------
        'INN_CLIP1' = 'integer',        # inner endpoint data
        'INN_POS1' = 'integer',
        'INN_SIDE1' = 'character',
        'INN_SEQ1' = 'character',
        #-------------
        'FLAG2' = 'integer',            # and again for the other side of the junction
        'RNAME2' = 'character',    
        'POS2' = 'integer',
        'MAPQ2' = 'integer',
        'SEQ2' = 'character',
        #-------------
        'OUT_CLIP2' = 'integer',        
        'OUT_POS2' = 'integer',
        'OUT_SIDE2' = 'character',
        'OUT_SEQ2' = 'character',
        #-------------
        'INN_CLIP2' = 'integer',            
        'INN_POS2' = 'integer',
        'INN_SIDE2' = 'character',
        'INN_SEQ2' = 'character',
        #-------------
        'MOL_OUT_POS1' = 'integer',    # outer position of molecule (not junction), only set for split junctions
        'MOL_OUT_POS2' = 'integer' 
    )
)
extract$clips <- c(
    extract$junctions[1:100 <= which(names(extract$junctions)=="OUT_SEQ1")],
    list(
        'FLAG2' = 'integer',     
        'RNAME2' = 'character',    
        'POS2' = 'integer',
        'MAPQ2' = 'integer',
        'SEQ2' = 'character',
        #-------------
        'OUT_CLIP2' = 'integer',        
        'OUT_POS2' = 'integer',
        'OUT_SIDE2' = 'character',
        'OUT_SEQ2' = 'character'
    )
)

#----------------------------------------------------------------------
# compile file columns
#----------------------------------------------------------------------
compile <- list(
    edges = list(
        'NODE1' = 'character',
        'NODE2' = 'character',
        'JXN_TYPE' = 'character',
        'TARGET_CLASS' = 'character',
        'MOL_IDS' = 'character',
        'COUNT' = 'integer',
        'COUNT_DISTINCT' = 'integer' 
    ),
    nodes = list(
        'NODE' = 'character', # node-level data
        'CLIP_LEN' = 'integer',
        'CLIP_SEQ' = 'character',
        'NODE_CLASS' = 'integer',
        #---------------
        'JXN_TYPE' = 'character', # edge/junction-level data
        'JXN_N' = 'integer',
        #---------------
        'FLAG' = 'integer', # alignment-level data
        'POS' = 'integer',
        'MAPQ' = 'integer',
        'CIGAR' = 'character',
        'SEQ' = 'character',
        'ALN_N' = 'integer',
        #---------------
        'UMI' = 'integer', # molecule-level data
        'MOL_ID' = 'integer',   
        'IS_DUPLEX' = 'integer',   
        'STRAND_COUNT1' = 'integer',
        'STRAND_COUNT2' = 'integer',
        'IS_MERGED' = 'integer',
        'TARGET_CLASS' = 'character', # based on molecule outer endpoints (not a junction)
        'MOL_CLASS' = 'character',
        'SHARED_PROPER' = 'integer',
        'IS_OUTER_CLIP1' = 'integer',
        'IS_OUTER_CLIP2' = 'integer',
        'IS_ORPHAN' = 'integer',
        'MOL_STRAND' = 'integer'
    )
)

#----------------------------------------------------------------------
# find, i.e. structural variant, columns
#----------------------------------------------------------------------
#find <- list(
#    molecules = c(list('SV_ID'='integer'), extract$junctions),
#    clips     = c(list('SV_ID'='integer'), extract$clips),
#    structural_variants = list(
#        'SV_ID' = 'numeric',
#        'MOL_IDS' = 'character', # comma-delimited
#        #-------------
#        'TARGET_CLASS' = 'character', # TT tt etc.
#        'JXN_TYPE' = 'character',  # DIPTL
#        'FLAG1' = 'integer',
#        'RNAME1' = 'character',
#        'INN_SIDE1' = 'character',
#        'FLAG2' = 'integer',
#        'RNAME2' = 'character',
#        'INN_SIDE2' = 'character',
#        #-------------
#        'PROX_OUT_POS1' = 'integer', # genome data 1
#        'PROX_JXN_POS1' = 'integer',
#        'PROX_JXN_POS_PLUS1' = 'integer',
#        'IS_CLIPPED1' = 'integer',
#        'JXN_JXN_POS1' = 'integer',
#        'JXN_MAPQ1' = 'integer',
#        #-------------
#        'PROX_OUT_POS2' = 'integer', # genome data 2
#        'PROX_JXN_POS2' = 'integer',
#        'PROX_JXN_POS_PLUS2' = 'integer',
#        'IS_CLIPPED2' = 'integer',
#        'JXN_JXN_POS2' = 'integer',
#        'JXN_MAPQ2' = 'integer',
#        #-------------
#        'CALL_LEN' = 'integer', # junction call data
#        'JXN_SEQ' = 'character',
#        'JXN_DEPTH' = 'character',
#        'MICROHOM_LEN' = 'integer',
#        'MICROHOM_MATCH' = 'integer',
#        'JXN_BASES' = 'character',
#        'SV_SIZE' = 'integer',
#        #-------------
#        'N_TOTAL' = 'integer', # junction evidence tallies
#        'N_SPLITS' = 'integer',
#        'N_GAPS' = 'integer',
#        'N_OUTER_CLIPS' = 'integer',
#        'N_DUPLEX' = 'integer',
#        'NET_STRAND_COUNT' = 'integer',
#        'N_SHARED_PROPER' = 'integer',
#        'N_UMI_PURGED' = 'integer'
#    ),
#    index = list( # binary index of files with supplemental information about SVs
#        'SV_ID' = 'numeric',
#        'MOL_OFFSET' = 'integer',
#        'MOL_LENGTH' = 'integer',
#        'CLIP_OFFSET' = 'integer',
#        'CLIP_LENGTH' = 'integer',   
#        'READS_OFFSET' = 'integer',
#        'READS_LENGTH' = 'integer',
#        'ALIGNS_OFFSET' = 'integer',
#        'ALIGNS_LENGTH' = 'integer',
#        'ALIGNS_PLUS_OFFSET' = 'integer',
#        'ALIGNS_PLUS_LENGTH' = 'integer'
#    )
#)
find <- list(
    structural_variants = list( 
        SV_ID = "integer", # SV identifiers
        JUNCTION_NAME = "character",  # the query junction
        MATCHING_NAMES = "character", # e.g. gaps that match a split query (should NOT be present as a JUNCTION_NAME in any row)
        OTHER_NAMES = "character",    # other non-matching SV junctions carried in the same molecules
        #-------------
        TARGET_CLASS = "character", # junction-level data (NB: this is a _new_ target class that applies to junction nodes
        JXN_TYPE = "character",
        #-------------
        CHROM_1 = "character", # node-level data
        SIDE_1 = "character",
        POS_1 = "integer",
        CHROM_2 = "character",
        SIDE_2 = "character",
        POS_2 = "integer",
        #-------------
        JXN_SEQ = "character", # joint (re)construction data
        MERGE_LEN = "integer",
        FAIDX_PADDING = "integer",
        GEN_REF_1 = "character",
        GEN_REF_2 = "character",
        #-------------
        MICROHOM_LEN = "integer", # joint data
        MICROHOM_MATCH = "integer",
        JXN_BASES = "character",
        SV_SIZE = "integer",     
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
        CHUNK_OFFSET = "integer", # index into the all_nodes evidence file
        CHUNK_SIZE = "integer"
    ),
    all_nodes = list(
        SV_ID = "integer", # numeric identifier for the SV junction
        NODE_N = "integer", # which SV junction node this matched
        IS_JUNCTION_NODE = "logical", # node is consistent with the query SV junction (false for non-SV nodes in molecules)
        IS_SEED_NODE = "logical", # node is an exact match to the the seed/query junction
        IS_REF_NODE = "logical", # this node was one of two used to characterize the junction
        IS_RC = "logical", # this SEQ was rc'ed as it was an inversion outer clip
        IS_REPEAT = "integer", # set in a later step to TRUE/1 if node was already claimed by a lower numbered SV
        N_COLLAPSED = "integer" # number of other initial molecules that were merged into this one as presumed duplicates
    )                      
)
find$all_nodes <- c(find$all_nodes, compile$nodes)

#----------------------------------------------------------------------
# compare, i.e. multi-sample structural variant, columns
#----------------------------------------------------------------------
compare <- list(
    structural_variants = c(list(
        'COMPARE_ID' = 'integer',
        'N_MATCHING_SVS' = 'integer',
        'N_MATCHING_SAMPLES' = 'integer',
        'HAS_COLLISION' = 'integer',
        'MATCHING_SAMPLES' = 'character',
        'SAMPLE' = 'character',
        'IS_COLLISION' = 'integer'
    ), find$structural_variants )
)

#----------------------------------------------------------------------
# molecule column output formatting
#----------------------------------------------------------------------
nodeTableNames <- list( # select and rename (for space) the output table columns
    MOL_ID='molId',JXN_N='jN',ALN_N='aN',
    NODE_N='nN',NODE_CLASS='nClss',JXN_TYPE='jTyp',NODE='node',CLIP_LEN='clip',CIGAR='cgr',MAPQ='mapQ',
    IS_JUNCTION_NODE='jNode',IS_REF_NODE='rNode',IS_SEED_NODE='sNode',IS_REPEAT='repeat',#UMI1='umi1',UMI2='umi2',
    IS_DUPLEX='dup',IS_MERGED='mrg',SHARED_PROPER='nPrp'
)
nodeToCommify <- c(
    #'innPos1','outPos1','innPos2','outPos2','molOutPos1','molOutPos2'
)
        
        
#----------------------------------------------------------------------
# structural variant column output formatting
#----------------------------------------------------------------------
#svTableNames <- list( # select and rename (for space) the output table columns
#    SAMPLE="smp",SV_ID='svId',
#    TARGET_CLASS='tgt',JXN_TYPE='SV',
#    #N_TOTAL='mol',
#    N_GAPS='gap',N_SPLITS='split',N_OUTER_CLIPS='clip',N_UMI_PURGED='purge',
#    N_DUPLEX='dup',NET_STRAND_COUNT='net',FRAC_SHARED_PROPER='fPrp',
#    RNAME1='chr1',PROX_JXN_POS1='pos1',JXN_MAPQ1='Q1', # ,INN_SIDE1='ori1'
#    RNAME2='chr2',PROX_JXN_POS2='pos2',JXN_MAPQ2='Q2', # ,INN_SIDE2='ori2'
#    SV_SIZE='size',
#    MICROHOM_LEN='uHom',JXN_BASES='jxnSeq'
#)
svTableNames <- list( # select and rename (for space) the output table columns
    SAMPLE="smp",SV_ID='svId',
    TARGET_CLASS='tgt',JXN_TYPE='SV',
    #N_TOTAL='mol',
    N_GAPS='gap',N_SPLITS='split',N_OUTER_CLIPS='clip',#N_UMI_PURGED='purge',
    N_DUPLEX='dup',NET_STRAND_COUNT='net',FRAC_SHARED_PROPER='fPrp',
    CHROM_1='chr1',POS_1='pos1',#JXN_MAPQ1='Q1', # ,INN_SIDE1='ori1'
    CHROM_2='chr2',POS_2='pos2',#JXN_MAPQ2='Q2', # ,INN_SIDE2='ori2'
    SV_SIZE='size',
    MICROHOM_LEN='uHom',JXN_BASES='jxnSeq'
)  
svToCommify <- c(
    'pos1','pos2','size'
)
svExcelColumns <- c(
    'SAMPLE', 'SV_ID',
    'TARGET_CLASS', 'JXN_TYPE', 'SV_SIZE', 'size',
    'N_SPLITS','N_GAPS','N_OUTER_CLIPS','N_DUPLEX','NET_STRAND_COUNT','FRAC_SHARED_PROPER',#,'N_UMI_PURGED'
    'RNAME1', 'PROX_JXN_POS1', 'INN_SIDE1','JXN_MAPQ1', 
    'RNAME2', 'PROX_JXN_POS2', 'INN_SIDE2','JXN_MAPQ2',
    'MICROHOM_LEN', 'MICROHOM_MATCH', 'JXN_BASES',   
    'PROX_OUT_POS1','JXN_JXN_POS1',
    'PROX_OUT_POS1','JXN_JXN_POS1',
    'JXN_SEQ', 'JXN_DEPTH'
)

#----------------------------------------------------------------------
# other lookups
#----------------------------------------------------------------------
jxnClasses <- list( # apply to molecules used as SV junction evidence
   'GAP'='0',
   'SPLIT'='1',
   'CLIP'='2'
)
jxnTypes = list( # possible chromosome orientations of the two sides of a junction
    D = "Dup",
    L = "Del",
    I = "Inv",
    T = "Trans",
    '?' = "Undet"
)
jxnTypes_peaks = list( # possible chromosome orientations of the two sides of a junction
    D = "Dup",
    L = "Del",
    I = "Inv"
)

