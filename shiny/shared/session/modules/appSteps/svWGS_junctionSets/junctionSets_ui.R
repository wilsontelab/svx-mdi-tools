#----------------------------------------------------------------------
# UI components for the junctionSets appStep module
#----------------------------------------------------------------------

# module ui function
svWGS_junctionSetsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$svWGS_junctionSets)

    # UI functions
    junctionAxisChoices <- list(
        "Junction CN"       = "junctionCN",
        "Outer Flank CN"    = "outerFlankCN",
        "Inner Flank CN"    = "innerFlankCN",
        "Event Span CN"     = "eventSpanCN",
        "Flank CNC"         = "flankCNC",
        "Event Span CNC"    = "eventCNC",
        "CNC"               = "chromCNC",
        "# of Samples"      = "N_SAMPLES",
        "log10(SV Size)"    = "SV_SIZE"
    )
    plotBox_ <- function(title, column1UI, column2UI, column3UI = NULL){
        box(
            title = title,
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            # collapsible = TRUE,
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
    jxnPlotUI <- function(id, title, xDefault, yDefault){
        plotBox_(
            title,
            tags$div(
                style = "max-height: 400px; overflow: auto;",
                radioButtons(ns(paste0(id, "YAxis")), "Y Axis", choices = junctionAxisChoices, selected = yDefault, inline = FALSE)
            ),
            interactiveScatterplotUI(ns(id), height = '400px'),
            tagList(
                jxnPlotBottomUI( radioButtons(ns(paste0(id, "XAxis")), "X Axis",     choices = junctionAxisChoices, selected = xDefault, inline = TRUE) )
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
        
        # junction exploration, with gating
        fluidRow(
            jxnPlotUI("jxnPlot",      "Junction Properties",       
                      xDefault = "outerFlankCN", yDefault = "flankCNC"),
            jxnPlotUI("gatedJxnPlot", "Junction Properties, Gated", 
                      xDefault = "junctionCN", yDefault = "SV_SIZE")
        ),

        # junction filtering
        fluidRow(
            box(
                width = 12,
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
            )
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
                "Gated SV Junctions",
                width = 12,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                downloadable = TRUE 
            )  
        ),

        # ... or expansion details about that object
        # row selection in expansionTable might populate expansionUI
        fluidRow(
            # class = "expansionTableWrapper browserExpansionWrapper",
            # style = "display: none;",
            bufferedTableUI(
                id = ns("expansionTable"), 
                title = "Supporting Molecules", 
                downloadable = TRUE, 
                width = 12,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE
            )  
        ),

        # a place for arbitrary, track-defined UI content based on expand2 or other click actions
        fluidRow(
            uiOutput(ns("expansionUI"))
        )
    )
}
