#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# objects instantiated here are available to all appStep modules in a session

# TODO: create mechanism for clearing cache via the UI
if(!exists('svxCache', envir = .GlobalEnv)){
    assign('svxCache', list(svs = list(), molecules = list()), envir = .GlobalEnv)
}
fillSvxCache <- function(sourceId, type, file){
    if(!is.null(svxCache[[type]][[file]])) return()
    svxCache[[type]][[file]] <<- {
        reportProgress("loading svx table")
        reportProgress(file)
        metadataFile <- getSourceFilePath(sourceId, "metadata")
        metadata <- read_yaml(metadataFile)
        metadata <- strsplit(unlist(metadata), "\\s+")
        fread(
            file,
            sep = "\t",        
            header = TRUE,
            stringsAsFactors = FALSE,
            colClasses = c(unname(unlist(SVX$find$structural_variants)), rep("integer", length(metadata$SAMPLES)))
        )
    }
}

# appServer function called after all modules are instantiated
appServer <- function(){


}
