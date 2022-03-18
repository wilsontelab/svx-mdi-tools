
# functions to help load and parse environment variables
# expects that calling script has executed:
#   env <- as.list(Sys.getenv())
checkEnvVar <- function(x){
    if(is.null(env[[x]])) stop(paste('missing environment variable:', x))
}
loadEnvVars <- function(vars = list(
    string  = character(), # var names by type
    integer = character(),
    double  = character()
)){
    invisible(sapply(unlist(vars), checkEnvVar))
    if(length(vars$integer) > 0) env[vars$integer] <- as.integer(env[vars$integer])
    if(length(vars$double)  > 0) env[vars$double]  <- as.double( env[vars$double])
}

# functions to help report summary values
reportStat <- function(value, message){
    message(paste(value, message, sep = "\t"))
}
