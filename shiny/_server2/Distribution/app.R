
#----------------------------------------------------------------------
# app.R assembles and launches the Shiny server application
#----------------------------------------------------------------------

# get environment variables
serverEnv <- as.list(Sys.getenv())

# source the app.R script common to several servers
source(paste(serverEnv$ACTIONS_PATH, '..', 'app.R', sep="/"))

