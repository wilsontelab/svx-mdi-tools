#----------------------------------------------------------------------
# server components for the svPoreSampleSelector widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svPoreSampleSelectorServer <- function(id, sampleSelection = "multiple") { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'svPoreSampleSelector'

#----------------------------------------------------------------------
# fill samples table when a source is selected
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single")
sourceSamplesReactive <- reactive({ 
    sourceId <- sourceId()
    req(sourceId)
    sourceSamplesFile <- expandSourceFilePath(sourceId, "sourceSamples.rds")
    if(!file.exists(sourceSamplesFile)){
        startSpinner(session, message = "reading source samples")
        x <- readRDS(getSourceFilePath(sourceId, "segmentsFile")) 
        x <- x[, .(nSegments = .N), keyby = .(sample)]
        saveRDS(x, sourceSamplesFile)
        stopSpinner(session)
    }
    readRDS(sourceSamplesFile)
})
samples <- bufferedTableServer(
    "samples",
    id,
    input,
    sourceSamplesReactive,
    selection = sampleSelection,
    selectionFn = function(selectedRows) NULL,
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
)

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
list(
    sourceId = sourceId, 
    samples = reactive({
        i <- samples$rows_selected()
        req(i)
        sourceSamplesReactive()[i, sample]
    })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
