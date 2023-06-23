#-------------------------------------------------------------------------------------
# svPore analyze constants
#-------------------------------------------------------------------------------------

# extract edges column definition
coreEdgeCols <- c( 
    "qName",
    "node1",
    "qStart",
    "node2",
    "qEnd",
    "mapQ",
    "cigar",
    "gapCompressedIdentity",
    "edgeType",
    "eventSize",
    "insertSize",
    "nStrands"
)
edgeExtensions1 <- c(
    "baseQual",
    "alnBaseQual",
    "alnSize",
    "sStart",
    "sEnd"    
)
edgeAdapterScores <- c(
    "clip5",
    "score5",
    "nBases5",
    "start5",
    "end5",
    "clip3",
    "score3",
    "nBases3",
    "start3",
    "end3"    
)
controlAdapterScores <- c(
    "score5C",
    "nBases5C",
    "start5C",
    "end5C",
    "score3C",
    "nBases3C",
    "start3C",
    "end3C"    
)
edgeExtensions2 <- c(
    "channel",
    "pod5File",
    "blockN",
    "edgeN",
    "readI"    
)
edgesCols <- c( 
    coreEdgeCols,
    edgeExtensions1,
    edgeAdapterScores,
    controlAdapterScores,
    edgeExtensions2
)
edgesColClasses <- c(
    "character",
    "integer64",
    "integer",
    "integer64",
    "integer",
    "integer",
    "character",
    "numeric",
    "character",
    "integer",
    "integer",
    "integer",
    #-------------,
    "numeric",
    "numeric",
    "integer",
    "integer",
    "integer",
    #-------------,
    "integer",
    "numeric",
    "integer",
    "integer",
    "integer",
    "integer",
    "numeric",
    "integer",
    "integer",
    "integer",
    "numeric",
    "integer",
    "integer",
    "integer",
    "numeric",
    "integer",
    "integer",
    "integer",
    #-------------,
    "integer",
    "character",
    "integer",
    "integer",
    "integer"
)

# edge types
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

# debug plot constants
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
