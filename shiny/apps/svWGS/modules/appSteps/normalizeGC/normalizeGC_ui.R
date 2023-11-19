#----------------------------------------------------------------------
# UI components for the normalizeGC appStep module
#----------------------------------------------------------------------

# module ui function
normalizeGCUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$normalizeGC)
    junctionAxisChoices <- list(
        "Junction CN"       = "junctionCN",
        "Non-SV Flank CN"   = "outerFlankCN",
        "SV Flank CN"       = "innerFlankCN",
        "SV Span CN"        = "eventSpanCN",
        "SV Flank CNC"      = "flankCNC",
        "SV Span CNC"       = "eventCNC",
        "CNC"               = "chromCNC",
        "Fraction GC"       = "gc",
        "# of Samples"      = "N_SAMPLES",              
        "log10(SV Size)"    = "SV_SIZE"
    )

    # UI functions
    plotBox_ <- function(title, column1UI, column2UI, column3UI = NULL){
        box(
            title = title,
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            column(
                width = 3,
                column1UI
            ),
            column(
                width = 9,
                column2UI
            ),
            column(
                width = 12,
                column3UI
            )
        )
    }
    jxnPlotBottomUI <- function(ui) tags$div(
        style = "margin-bottom: 10px;",
        ui
    )
    jxnPlotUI <- function(id, title, xDefault, yDefault, filters = FALSE){
        plotBox_(
            title,
            tags$div(
                style = "max-height: 400px; overflow: auto;",
                radioButtons(ns(paste0(id, "YAxis")), "Y Axis", choices = junctionAxisChoices, selected = yDefault, inline = FALSE)
            ),
            interactiveScatterplotUI(ns(id), height = '400px'),
            tagList(
                jxnPlotBottomUI( radioButtons(ns(paste0(id, "XAxis")), "X Axis",     choices = junctionAxisChoices, selected = xDefault, inline = TRUE) ),
                if(filters) fluidRow(
                    column(
                        width = 3,
                        jxnPlotBottomUI( radioButtons(ns("jxnPlotNSamples"), "# of Samples", 
                                                    choices = c("Unique", "Shared"), selected = "Unique", inline = TRUE) )                        
                    ),
                    column(
                        width = 3,
                        jxnPlotBottomUI( checkboxGroupInput(ns("jxnPlotSVTypes"), "SV Types", 
                                                            choices  = c("Del","Dup","Inv","Trans"), 
                                                            selected = c("Del","Dup"), inline = TRUE))                  
                    ),
                    column(
                        width = 3,
                        jxnPlotBottomUI( numericInput(ns("jxnPlotMinNInstances"), "Min # of Reads", value = 3, min = 1, max = 100, step = 1))           
                    ),
                    column(
                        width = 3,
                        jxnPlotBottomUI( numericInput(ns("jxnPlotMinSequenced"), "Min Sequenced Reads", value = 2, min = 0, max = 100, step = 1))           
                    )
                ) else NULL
            )
        )
    }

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

        # saved junction sets (save button is at the bottom...)
        fluidRow(
            bufferedTableUI(
                ns("savedJunctionSets"), 
                title = "Saved Junction Sets", 
                width = 12,
                collapsible = TRUE,
                collapsed = FALSE,
                downloadable = TRUE 
            )
        ),

        # data source selectors
        fluidRow(
            dataSourceTableUI(
                ns("source"), 
                "Data Source", 
                width = 6, 
                collapsible = FALSE,
                inFluidRow = FALSE
            ),
            bufferedTableUI(
                ns("sampleGenome"),
                "Sample with Genome",
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = FALSE
            )            
        ),

        # GC fit and chromosome profile plots
        fluidRow(
            plotBox_(
                "GC Bias Fit",
                tags$div(
                    numericInput(ns("ploidy"),     "Median Ploidy", value = 2, min = 1, max = 4, step = 1),
                    selectInput(ns("expectedSex"), "Expected Sex",  selected = "auto", choices = c("auto","XX","XY"))
                ),
                interactiveScatterplotUI(ns("gcBiasPlot"), height = '400px')
            ),
            plotBox_(
                title = "Chromosome Profiles",
                tags$div(
                    style = "max-height: 400px; overflow: auto;",
                    radioButtons(ns("densityChrom"), "Chromosome", choices = c("genome"), inline = FALSE)
                ),
                interactiveScatterplotUI(ns("chromDensityPlot"), height = '400px')
            )
        ),
    
        # junction exploration, with gating
        fluidRow(
            jxnPlotUI("jxnPlot",      "Junction Properties",       
                      xDefault = "outerFlankCN", yDefault = "flankCNC", filters = TRUE),
            jxnPlotUI("gatedJxnPlot", "Junction Properties, Gated", 
                      xDefault = "junctionCN", yDefault = "flankCNC")
        ),

        # saving of gated SV junction sets (list of saved sets at the top...)
        fluidRow( 
            style = "margin-bottom: 10px;",
            column(
                offset =4,
                width = 4,
                bsButton(
                    ns("saveJunctionSet"),
                    "Save Junction Set",
                    style = "success",
                    block = TRUE
                )
            ),
            column(
                width = 4,
                uiOutput(ns('saveJunctionSetFeedback'))
            )
        ),

        # table of gated junctions as a candidate junction set
        fluidRow(
            bufferedTableUI(
                ns("gatedSvs"),
                "Gated SVs",
                width = 12,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                downloadable = TRUE 
            )  
        ),

        NULL
    )
}
