#-------------------------------------------------------------------------------------
# svPore analyze constants
#-------------------------------------------------------------------------------------
edgesCols <- c( 
    "qName",
    "node1",
    "cigar", # omit "1" suffix since will purge the redudancy to a single entry per alignment
    "blastIdentity",
    "gapCompressedIdentity",
    "node2",
    "cigar2",
    "blastIdentity2",
    "gapCompressedIdentity2",
    "edgeType",
    "mapQ",
    "eventSize",
    "insertSize",
    "xStart",
    "xEnd",
    "edgeClass",
    "nStrands"
)
edgesColClasses <- c(
    "character",
    "integer",
    "character",
    "numeric",
    "numeric",
    "integer",
    "character",
    "numeric",
    "numeric",
    "character",
    "integer",
    "integer",
    "integer",
    "integer",
    "integer",
    "integer",
    "integer"
)
readsCols = c(
    "molType", 
    "qName",
    "qSeq"
    # ,
    # "qQual"
)
readsColClasses = c(
    "character",
    "character",
    "character"
    # ,
    # "character"
)

edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    PROPER        = "P"
)
inlineSvTypes <- c(edgeTypes$DELETION, edgeTypes$INSERTION)
incrementRefTypes <- c(edgeTypes$ALIGNMENT, edgeTypes$DELETION)
splitTypes <- c(edgeTypes$TRANSLOCATION, edgeTypes$INVERSION)

edgeClasses <- list(
    FROM_SPLIT = 2, # thus, the value implies the number of alignments 
    FROM_CIGAR = 1
)

plotlyColors <- list(
    blue    = '#1f77b4',  # muted blue
    orange  = '#ff7f0e',  # safety orange
    green   = '#2ca02c',  # cooked asparagus green
    red     = '#d62728',  # brick red
    purple  = '#9467bd',  # muted purple
    brown   = '#8c564b',  # chestnut brown
    pink    = '#e377c2',  # raspberry yogurt pink
    gray    = '#7f7f7f',  # middle gray
    yellow  = '#bcbd22',  # curry yellow-green
    teal    = '#17becf',  # blue-teal
    black   = 'black',
    grey    = '#7f7f7f'
)
edgeTypeColors <- list(
    "A" = rgb(0, 0, 0, 0.1), # the single type for a contiguous aligned segment

    "T" = rgb(1, 0, 0, 0.1), # edge/junction types (might be several per source molecule)
    "V" = rgb(0, 0.8, 0, 0.1),
    "U" = rgb(0.8, 0.8, 0, 0.1),
    "D" = rgb(0, 0, 1, 0.1),

    "?" = rgb(0, 0, 0, 0.1),
    "I" = rgb(0, 0, 0, 0.1), 
    "M" = rgb(0, 0, 0, 0.1),
    "P" = rgb(0, 0, 0, 0.1),
    "F" = rgb(0, 0, 0, 0.1),
    "R" = rgb(0, 0, 0, 0.1),
    "Q" = rgb(0, 0, 0, 0.1)
)
