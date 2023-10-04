#----------------------------------------------------------------------
# appServer() is called in session context and thus has access to:
#   input, output, session objects
#   values returned from app step modules
#----------------------------------------------------------------------

# genome browser support
library(bit64)
sessionCache <- new_dataCache('svAmpliconCache')

# constants
CONSTANTS$edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    MERGE_FAILURE = "M",
    PROPER        = "P"
)
CONSTANTS$junctionTypes <- list(
    D   = "deletion (D)",
    I   = "insertion (I)",
    U   = "duplication (U)",    
    V   = "inversion (V)",
    T   = "translocation (T)"
)
CONSTANTS$groupedPathClasses <- list(
    "-"  = "no SV",
    D   = "deletion (D)",
    I   = "insertion (I)",
    U   = "duplication (U)",    
    "V:V"  = "full inversion (VV)",
    V   = "partial inversion (V)",
    "T:T"  = "full translocation (TT)",
    T   = "partial translocation (T)",
    "D:I"  = "paired indel (DI,ID)",
    "I:D"  = "paired indel (DI,ID)"
    # all other path classes listed as "other"
)

# appServer() is called after all modules are instantiated
appServer <- function(){

    # objects instantiated here are available to this app step only

}
