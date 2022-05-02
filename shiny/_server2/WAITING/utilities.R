
#----------------------------------------------------------------------
# utilities.R has universal server and R helper functions
#----------------------------------------------------------------------

# code development utilities for debugging and finding slow steps
verbose <- FALSE
reportProgress <- function(message){
    if(verbose) message(message)
}

# miscellaneous
`%notin%` <- Negate(`%in%`)

