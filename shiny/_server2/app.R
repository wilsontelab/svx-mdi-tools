
#----------------------------------------------------------------------
# app.R assembles and launches the Shiny server application
#----------------------------------------------------------------------

# load server-specific dependencies not found in _shiny/app.R
# none

# define the server initialization sequence
initalizeApp <- function(){
    
    # STEP 1: source server-specific scripts, typically needed for STEP 2
    # these are server level scripts scoped for all users (session scripts are sourced in server.R)
    sourceScripts(c(
        "constants.R",  
        "column_definitions.R",
        "data_sources_server.R"
    ), subDir="..") # ../<scipt> pattern >> shared by multiple server apps used by a pipeline
    
    # STEP 2: initialize variables required to configure the page
    setProjects()
    setSamples()
    
    # STEP 3: source the shiny server script that construct the page
    sourceScripts(c(
        "ui.R",
        "server.R"  
    ))
}

# source _shiny/app.R to load and launch server
source(paste(serverEnv$SHINY_PATH, 'app.R', sep="/"))

