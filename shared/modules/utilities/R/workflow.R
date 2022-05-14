# load, parse and save environment variables
# expects that calling script has executed:
#   env <- as.list(Sys.getenv())
checkEnvVar <- function(x){
    if(is.null(env[[x]])) stop(paste('missing environment variable:', x))
}
checkEnvVars <- function(vars = list(
    string  = character(), # var names by type
    integer = character(),
    double  = character(),
    logical = character()
)){
    invisible(sapply(unlist(vars), checkEnvVar))
    if(length(vars$integer) > 0) env[vars$integer] <<- as.integer(env[vars$integer])
    if(length(vars$double)  > 0) env[vars$double]  <<- as.double( env[vars$double])
    if(length(vars$logical) > 0) env[vars$logical] <<- as.logical(as.integer(env[vars$logical]))
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

# report and retrieve summary values
reportStat <- function(value, message){
    message(paste(value, message, sep = "\t"))
}
getStatValue <- function(projectDir, sample, action, scriptId, statName, pipeline = NULL){
    if(is.null(pipeline)) pipeline <- env$PIPELINE_NAME
    dir <- file.path(projectDir, sample, pipeline, action, 'logs')
    if(!dir.exists(dir)) return(NA)
    fileName <- paste(sample, scriptId, "log", "txt", sep = ".")
    statsFile <- file.path(dir, fileName)
    if(!file.exists(statsFile)) return(NA)
    stats <- fread(statsFile, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if(nrow(stats) == 0) return(NA)
    names(stats) <- c('scriptId_', 'value_', 'statName_', 'description_')
    row <- stats[, last(.I[statName_ == statName])]
    if(length(row) == 0) NA else stats[row, value_]
}
