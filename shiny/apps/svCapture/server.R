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
    assign('svxCache', list(
        svs       = list(), 
        molecules = list(),
        targets   = list()
    ), envir = .GlobalEnv)
}
fillSvxCache <- function(sourceId, type, contentFileType, colClasses){
    file <- getSourceFilePath(sourceId, contentFileType)
    if(!is.null(svxCache[[type]][[file]])) return(file)
    svxCache[[type]][[file]] <<- {
        reportProgress("loading svx data")
        reportProgress(file)
        fread(
            file,
            sep = "\t",        
            header = TRUE,
            stringsAsFactors = FALSE,
            colClasses = colClasses()
        )
    }
    file
}
fillSvxCache_svs <- function(sourceId){
    fillSvxCache(sourceId, "svs", "structuralVariants", function(){
        metadataFile <- getSourceFilePath(sourceId, "metadata")
        metadata <- read_yaml(metadataFile)
        metadata <- strsplit(unlist(metadata), "\\s+") 
        c(unname(unlist(SVX$find$structural_variants)), rep("integer", length(metadata$SAMPLES)))       
    })
}
fillSvxCache_targets <- function(sourceId){
    fillSvxCache(sourceId, "targets", "targetsBed", function() NULL)
}

# appServer function called after all modules are instantiated
appServer <- function(){


}
