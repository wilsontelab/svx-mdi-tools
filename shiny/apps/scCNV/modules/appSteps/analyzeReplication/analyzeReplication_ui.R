#----------------------------------------------------------------------
# UI components for the analyzeReplication appStep module
#----------------------------------------------------------------------

# module ui function
analyzeReplicationUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$analyzeReplication)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        # settings = TRUE,

        # appStep UI elements, populate as needed
        "module contents pending"
    )
}
