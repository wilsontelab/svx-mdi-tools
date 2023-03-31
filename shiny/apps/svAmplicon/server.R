#----------------------------------------------------------------------
# appServer() is called in session context and thus has access to:
#   input, output, session objects
#   values returned from app step modules
#----------------------------------------------------------------------

# objects instantiated here are available to all appStep modules in a session
CONSTANTS$edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    MERGE_FAILURE = "M",
    PROPER        = "P",
    FUSED_MERGE_FAILURE = "F",
    REJECTED_INDEL = "R",
    FUSED_MERGE_FAILURE_REJECTED_INDEL = "Q"
)

# appServer() is called after all modules are instantiated
appServer <- function(){

    # objects instantiated here are available to this app step only

}
