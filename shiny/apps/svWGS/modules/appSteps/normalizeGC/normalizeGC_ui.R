#----------------------------------------------------------------------
# UI components for the normalizeGC appStep module
#----------------------------------------------------------------------

# module ui function
normalizeGCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalizeGC)

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

        # data source selector
        dataSourceTableUI(
            ns("source"), 
            "Sample Source", 
            width = 12, 
            collapsible = FALSE
        ),

        # one row per source sample for adjusting GC bias normalization
        uiOutput(ns("samples"))
    )
}
