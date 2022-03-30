
#----------------------------------------------------------------------
# ui.R defines the html page layout and defines user inputs
#----------------------------------------------------------------------

# number of samples that can be concurrently selected
nSamples <- 5 

# page layout
ui <- fluidPage(
    
    # some page styles
    style="padding: 10px; font-size: 12px;",
    tags$head(tags$style(HTML( # hack fixes stupid default Shiny checkbox group indenting
        ".checkbox-inline, .radio-inline  { 
            margin-left: 4px;
            margin-right: 4px;
        }
        .checkbox-inline+.checkbox-inline, .radio-inline+.radio-inline {
            margin-left: 4px;
            margin-right: 4px;
        }
        hr {
            margin-top: 0;
            margin-bottom: 10px;
            border: 0;
            border-top: 1px solid #666;
        }
        #downloadExcel {
            color: rgb(0,0,200);
        }
        "
    ))),
    
    # standard Shiny main page layout
    sidebarLayout(
        
        # top level inputs that control the output
        sidebarPanel(

            # sample selection
            selectInput(inputId = 'project',
                        label = 'Project',
                        choices = projects),
            
            # filters for selecting the plotted and tabulated SVs
            selectInput(inputId = 'captureTarget',
                        label = 'Capture Target',
                        choices = c('-')),
            radioButtons(inputId = "duplexFilter",
                         label = "Has Duplex Molecule",
                         inline = TRUE,
                         choices = c('all','yes','no')),
            radioButtons(inputId = "netReadPairFilter",
                         label = "Net Read Pairs",
                         inline = TRUE,
                         choices = c('all','1','>3','>10')),
            radioButtons(inputId = "fracProperEnds",
                         label = "Avg. Shared w/Proper",
                         inline = TRUE,
                         choices = c('all','<1','>=1')),            
            radioButtons(inputId = "moleculeCountFilter",
                         label = "Uniq. Molecule Count",
                         inline = TRUE,
                         choices = c('all','1','>1','>3')),
            radioButtons(inputId = "sampleCountFilter",
                         label = "Matching Sample Count",
                         inline = TRUE,
                         choices = c('all','1','>1')),
            checkboxGroupInput(inputId = 'targetTypeFilter',
                               label = "Target End Pairing", inline="TRUE",
                               choices  = targetPairTypes,
                               selected = defaultTargetPairTypes),
            radioButtons(inputId = "svTypeFilter",
                         label = "SV Type",
                         inline = TRUE,
                         choiceNames  = c('all', unname(unlist(jxnTypes))),
                         choiceValues = c('all', unname(names(jxnTypes))) ),
            textInput(inputId = "minSvSizeFilter",
                         label = "Min. SV Size",
                         value  = 0 ),
            radioButtons(inputId = "sequencedFilter",
                         label = "Junction Sequenced",
                         inline = TRUE,
                         choices = c('all','yes','no')),
            
            ## some preset for common filter combinations
            #actionButton(inputId = "candidateSV", "candSV", inline=TRUE),
            #actionButton(inputId = "ligationArtifact", "ligArt", inline=TRUE),
            #actionButton(inputId = "chimericPCR", "chPCR", inline=TRUE),
            #br(), br(),
            
            # parameters influencing the appearance of the plots
            hr(),
            checkboxGroupInput(inputId = 'showPlots',
                               label = "Show Plots", inline="TRUE",
                               choices  = c('triangle','linear','corr','size_uHom','uHom'),
                               selected = c('triangle','size_uHom')),
            radioButtons(inputId = "plotColors",
                         label = "Plot Color Scheme",
                         inline = TRUE,
                         choices  = names(plotColors)),            
            radioButtons(inputId = "pointCex",
                         label = "Point Size (R cex)",
                         inline = TRUE,
                         choices = c('1','0.8','0.5')),

            # parameters influencing the appearance of the sequence alignments
            hr(),
            radioButtons(inputId = "clipMode",
                         label = "Color Clips",
                         inline = TRUE,
                         choices  = c('base','grey','none')),
            #radioButtons(inputId = "jxnDisplayType",
            #             label = "Junction Display",
            #             inline = TRUE,
            #             choices  = c('Genome','Genome+','Reads')),
            #radioButtons(inputId = "charPerLine",
            #             label = "Bases Per Line",
            #             inline = TRUE,
            #             choices  = seq(100,200,50)),
            
            # some action buttons
            hr(),
            actionButton(inputId = "clearSamples", "Clear Samples", inline=TRUE),
            actionButton(inputId = "resetFilters", "Reset Filters"),
            #actionButton(inputId = "refreshServer", "Refresh Server"),

            # width of the option section
            width=3
        ),

        # the output elements top to bottom
        mainPanel(
            width=9,
            fluidRow(
                column(1, "Sample(s)", style="font-weight: bold;"),
                lapply(1:nSamples, function(i){
                    column(2, selectInput(inputId = paste0('sample_', i),
                        label = NULL, choices = samples))                   
                }),
                column(1, downloadButton('downloadExcel', label="Excel"))
            ),                
            plotOutput(outputId="plotOutput"), #, height='900px'
            withSpinner(dataTableOutput(outputId="svTableOutput")),
            tags$div(
                id="junctionPlotDiv",
                imageOutput("consensusPlot", width="100%", height="auto"),
                br(),
                imageOutput("junctionPlot",  width="100%", height="auto"),
                width="100%",
                height="auto",
                style="overflow-x: scroll;"
            ),
            #uiOutput("junctionPlot")
            #plotOutput(outputId="zoomPlots", height='250px'),            
            #span(
            #     htmlOutput(outputId="junctionZoom"),
            #     style="font-family: monospace; font-weight: bold;"
            #),
            dataTableOutput(outputId="nodeTableOutput"),
            ""
        )
    )
)

