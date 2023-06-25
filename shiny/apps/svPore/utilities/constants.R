#-------------------------------------------------------------------------------------
# svPore constants
#-------------------------------------------------------------------------------------
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
junctionTypes <- list(
    All           = "X",
    Translocation = "T", # edge/junction types (might be several per source molecule)
    Inversion     = "V",
    Duplication   = "U",
    Deletion      = "D"
)
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
junctionTypeLines <- list(
    "T" = 4,
    "V" = 3,
    "U" = 2,
    "D" = 1,
    "I" = 0,    
    "?" = -1
)
