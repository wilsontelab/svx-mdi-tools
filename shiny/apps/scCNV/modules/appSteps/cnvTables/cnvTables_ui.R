#----------------------------------------------------------------------
# UI components for the cnvTables appStep module
#----------------------------------------------------------------------

# module ui function
cnvTablesUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$cnvTables)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        collapsibleDivUI(
            ns("leaderText"), 
            "Select one or more project sources to construct tables that summarize your kept and grouped CNVs across all selected samples.",
            "Successive tables are constructed with increasing levels of grouping. The first table has one row per CNV per cell, the second",
            "has one row per grouped CNV with sample and cell counts, the last table has cell counts per sample for each grouped CNV.",
            "In the last tables, the cell with the highest resolution is used to set the CNV boundaries."                              
        ),

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = TRUE,

        # box for selecting sample source
        scCnvDataSourceTableUI(ns),

        # sample selector
        sampleSelectorDivUI(ns),

        fluidRow(
            bufferedTableUI(
                ns("cnvsTable"), 
                title = "CNVs Table (one CNV+cell per row)", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = TRUE
            )
        ),

        fluidRow(
            bufferedTableUI(
                ns("groupsTable"), 
                title = "CNV Groups Table (aggregate cells per row)", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),

        fluidRow(
            bufferedTableUI(
                ns("groupsPivotTable"), 
                title = "CNV Groups Pivot (one sample per column)", 
                downloadable = TRUE,
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),
        NULL
    )
}
