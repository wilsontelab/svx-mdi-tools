#----------------------------------------------------------------------
# reactive components for the DESeq2 analysis type results viewer
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
DESeq2Server <- function(id, parentId, options, job, settings) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        parentNS <- function(id) paste(parentId, ns(id), sep = "-")
        module <- 'DESeq2' # for reportProgress tracing   
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# convenience accessors to settings and data 
#----------------------------------------------------------------------
plotThreshold <- function(name) settings$get('Plot_Thresholds', name)
sampleData <- reactive({ # the raw counts input; a matrix of [features,samples]
    job <- job()
    req(job)
    job$output$results$input$colData
})
counts <- reactive({ # the raw counts input; a matrix of [features,samples]
    job <- job()
    req(job)
    job$output$results$input$countData
})
groups <- reactive({
    job <- job()
    req(job)
    x <- getSampleSetAssignments(job$schema$Sample_Set, categoryNames = TRUE)
    categoryField <- job$schema$Analyze_By
    categoryNameField <- paste0(categoryField, 'Name')
    x$category <- x[[categoryField]]
    x$categoryName <- x[[categoryNameField]]
    x[, sampleId := paste(Project, Sample_ID)]
    x[, replicate := as.integer(factor(sampleId)), by = .(category)]
    x[, .(category, categoryName, replicate)]
})
de <- reactive({ # DESeq2 output; see str at bottom
    job <- job()
    req(job)
    job$output$results$output
})
significant <- reactive({ # logical vector for adjusted significance of each gene per settings threshold
    de <- de()
    req(de)
    ifelse(is.na(de$padj), FALSE, -log10(de$padj) >= plotThreshold('Log10_P_Adjusted') )
})

#----------------------------------------------------------------------
# user outcome reactives
#----------------------------------------------------------------------
marked <- reactiveVal(FALSE)

#----------------------------------------------------------------------
# gene results table
#----------------------------------------------------------------------
deCols <- list(
    row = 'Gene/Feature',
    baseMean = 'Mean_Count',
    log2FoldChange = 'Log2_Fold_Change',
    lfcSE = 'Std_Error',
    stat = 'Statistic',
    pvalue = 'P-Value',
    padj = 'P-Adjusted' 
)
numericCols <- names(deCols)[names(deCols) != 'row']
genesTableData <- reactive({
    de <- de()
    req(de)
    for (col in numericCols) de[[col]] <- signif(de[[col]], digits = 3)
    colnames(de) <- unlist(deCols[colnames(de)])
    marked( rep(FALSE, nrow(de)) )
    de    
})
genesTable <- resultsTableServer(
    'genesTable',
    parentId = id,
    type = 'long',
    tableData = genesTableData,
    getFileName = function(){ # no extension, module adds '.csv' and replaces spaces
        paste(getSchemaName(job()$schemaId), 'byGene', sep = '.')
    },
    marks = marked,
    parentNS = parentNS
)

#----------------------------------------------------------------------
# composite plots of all genes
#----------------------------------------------------------------------
getGenePlotData <- reactive({
    de <- de()
    req(de)
    x <- data.table(
        log2FoldChange = de$log2FoldChange,
        pvalue = -log10(de$pvalue),
        baseMean = log10(de$baseMean),
        gene = de$row,
        labels = ifelse(marked(), de$row, NA)
    )
    selected <- genesTable$selected()
    if(is.na(selected)) selected <- 0
    x[, trace := ifelse(.I == selected, 4, ifelse(significant(), 2, 1)) - 1] # use plotly standard trace colors
    x   
})
plotGeneClick <- function(gene){
    observeEvent(gene$clicked(), {
        i <- which(de()$row == gene$clicked()$customdata)
        genesTable$selectRow( i )
    })  
}
interactiveGenePlot <- function(id, xcol, ycol, xtitle, ytitle, labelDirs=list(x = 1, y = 1)){
    interactiveScatterplotServer(
        id,
        reactive({
            genePlotData <- copy( getGenePlotData() )
            req(genePlotData)
            genePlotData$x <- genePlotData[[xcol]]
            genePlotData$y <- genePlotData[[ycol]]
            genePlotData
        }), 
        accelerate = TRUE,
        xtitle = xtitle,
        ytitle = ytitle,
        clickable  = TRUE, 
        overplot = 'trace',
        overplotPointSize = 5,
        hoverText = 'gene',
        keyColumn = 'gene',
        labelCol = 'labels',
        labelDirs = labelDirs
    ) %>% plotGeneClick()
}
#----------------------------------------------------------------------
interactiveGenePlot(
    id = 'volcanoPlot',
    xcol = "log2FoldChange",
    ycol = "pvalue",
    xtitle = "Log2 Fold Change",
    ytitle = "Log10 P Value",
    labelDirs = list(x = -1, y = 1)
)
interactiveGenePlot(
    id = 'MAPlot',
    xcol = "baseMean",
    ycol = "log2FoldChange",
    xtitle = "Log10 Base Mean",
    ytitle = "Log2 Fold Change"
)

#----------------------------------------------------------------------
# gene results table and single-gene zoom from row click
#----------------------------------------------------------------------
selectedGeneName <- reactive({
    selectedRow <- genesTable$selected()
    req(selectedRow)    
    geneName <- de()[selectedRow, 'row']    
})
genePlotData <- reactive({ # temporary crude, much more can be done
    counts <- counts()
    req(counts)
    row <- which( rownames(counts) == selectedGeneName() )    
    dt <- cbind(data.table(value = counts[row, ]), groups())[order(-category)]
    dt[, ':='(
        group = categoryName,
        subgroup = replicate
    )]
    dt
})
interactiveBarplotServer(
    'genePlot',
    genePlotData,
    orientation = 'horizontal',
    xtitle = reactive({ paste(selectedGeneName(), "Raw Count") })
)

#----------------------------------------------------------------------
# set return value to job, to be used in developerTools
#----------------------------------------------------------------------
NULL

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

 #$ schema  :List of 6
 # ..$ Sample_Set   : chr "602a7b4d09b1484f72478a5badab3c74"
 # ..$ Analysis_Type: chr "DESeq2"
 # ..$ Analyze_By   : chr "Category1"
 # ..$ optionsHtml  : chr "Analyze_By = Group"
 # ..$ name         : chr "Analysis #1"
 # ..$ status       : int 0
 #$ schemaId: chr "8aa1fd62ff1c1261ba1f7f234e14dea3"
 #$ results :List of 3
 # ..$ status : num 0
 # ..$ message: chr "OK"
 # ..$ results:List of 2
 # .. ..$ input :List of 2
 # .. .. ..$ colData  :'data.frame':     4 obs. of  5 variables:
 # .. .. .. ..$ Source_ID: chr [1:4] "a2903be5721de792d7670b94efae0963" ...
 # .. .. .. ..$ Project  : chr [1:4] "7-LZ" "7-LZ" "7-LZ" "7-LZ"
 # .. .. .. ..$ Sample_ID: chr [1:4] "7-LZ-1" "7-LZ-2" "7-LZ-3" "7-LZ-4"
 # .. .. .. ..$ Category1: Factor w/ 2 levels "1","2": 1 2 1 2
 # .. .. .. ..$ Category2: Factor w/ 1 level "1": 1 1 1 1
 # .. .. ..$ countData: num [1:24062, 1:4] 0 497 74 4 113 0 480 277 61 205 ...
 # .. .. .. ..- attr(*, "dimnames")=List of 2
 # .. .. .. .. ..$ : chr [1:24062] "0610005C13Rik" "0610007P14Rik" "0610009B22Rik" "0610009L18Rik" ...
 # .. .. .. .. ..$ : chr [1:4] "7-LZ:7-LZ-1" "7-LZ:7-LZ-2" "7-LZ:7-LZ-3" "7-LZ:7-LZ-4"
 # .. ..$ output:'data.frame':   24062 obs. of  7 variables:
 # .. .. ..$ row           : chr [1:24062] "0610005C13Rik" "0610007P14Rik" "0610009B22Rik" "0610009L18Rik" ...
 # .. .. ..$ baseMean      : num [1:24062] 0 591.53 58.92 6.15 136.48 ...
 # .. .. ..$ log2FoldChange: num [1:24062] NA 0.0576 -0.8832 0.5079 0.4222 ...
 # .. .. ..$ lfcSE         : num [1:24062] NA 0.317 0.508 1.356 0.404 ...
 # .. .. ..$ stat          : num [1:24062] NA 0.181 -1.739 0.375 1.045 ...
 # .. .. ..$ pvalue        : num [1:24062] NA 0.8561 0.0821 0.708 0.296 ...
 # .. .. ..$ padj          : num [1:24062] NA 1 1 1 1 ...
  