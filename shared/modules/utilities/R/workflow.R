
# load, parse and save environment variables
# expects that calling script has executed:
#   env <- as.list(Sys.getenv())
checkEnvVar <- function(x){
    if(is.null(env[[x]])) stop(paste('missing environment variable:', x))
}
checkEnvVars <- function(vars = list(
    string  = character(), # var names by type
    integer = character(),
    double  = character()
)){
    invisible(sapply(unlist(vars), checkEnvVar))
    if(length(vars$integer) > 0) env[vars$integer] <<- as.integer(env[vars$integer])
    if(length(vars$double)  > 0) env[vars$double]  <<- as.double( env[vars$double])
}
writeEnvJson <- function(prefix){
    jsonlite::write_json(env, paste(prefix, 'environment', 'json', sep = "."))
}

# source scripts into the environment of the parent of this function (typically .GlobalEnv)
sourceScripts <- function(dir, scripts){
    suffix <- '.R'
    local <- parent.env(environment())
    for(script in scripts){
        if(!endsWith(script, suffix)) script <- paste0(script, suffix)
        source(file.path(dir, script), local = local)
    }
}

# report summary values
reportStat <- function(value, message){
    message(paste(value, message, sep = "\t"))
}
