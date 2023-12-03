#----------------------------------------------------------------------
# server components for svWGS explorer data selectors
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svWGS_explorer_dataSelectorsServer <- function(
    id,
    stepId,
    stepSettings, 
    sourceSelection = "multiple",
    sampleGenomesSelection = "multiple"
) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# selected data sources and associated source-level data
#----------------------------------------------------------------------
sourceIds <- dataSourceTableServer("sources", selection = sourceSelection) 
samples <- reactive({ # vector of the names of all co-analyzed samples
    sourceIds <- sourceIds()
    req(sourceIds)
    do.call(rbind, lapply(sourceIds, function(sourceId){
        source <- getSourceFromId(sourceId)
        source$manifest[, .(sourceId = sourceId, Project = Project, Sample_ID = Sample_ID)]
    }))
})
genomes <- reactive({ # data.table
    sourceIds <- sourceIds()
    req(sourceIds)
    do.call(rbind, lapply(sourceIds, function(sourceId){
        source <- getSourceFromId(sourceId)
        sourceGenome <- source$config$task$find$genome$genome
        customGenome <- loadCustomGenomeMetadata(sourceGenome)
        isCompositeGenome <- !is.null(customGenome) && isCompositeGenome2(customGenome)
        compositeDelimiter <- if(isCompositeGenome) getCustomCompositeDelimiter(customGenome) else NA
        genomes <- readRDS(getSourceFilePath(sourceId, "genomesFile"))
        genomes[, ":="(
            sourceId = sourceId,
            sourceGenome = sourceGenome,
            isCompositeGenome = isCompositeGenome,
            compositeDelimiter = compositeDelimiter
        )]
    }))
})
getJunctionGenome <- function(genome, sourceId_, chrom){
    if(genome$isCompositeGenome) sapply(chrom, function(x) strsplit(x, genome$compositeDelimiter)[[1]][2]) 
    else genome$sourceGenome
}
junctions <- reactive({
    sourceIds <- sourceIds()
    req(sourceIds) 
    samples <- samples()
    genomes <- genomes() # TODO: cache this?
    startSpinner(session, message = "loading junctions")
    x <- do.call(rbind, lapply(sourceIds, function(sourceId_){
        samples_ <- samples[sourceId == sourceId_]
        genome <- genomes[sourceId == sourceId_][1]
        x <- svx_loadJunctions(sourceId_, svWGS_loadJunctions)[, ":="(
            X = 1,
            sourceId = sourceId_,
            sourceGenome = genome$sourceGenome,
            Project = samples_$Project[1],
            Sample_ID = trimws(gsub(",", " ", samples)),
            Clonal = ifelse(nInstances > 1, "Clonal", "Single"),
            SV_Type = svx_jxnType_codeToX(edgeType, "name"),
            genome1 = getJunctionGenome(genome, sourceId_, cChrom1),
            genome2 = getJunctionGenome(genome, sourceId_, cChrom2)
        )][, ":="(
            Genome = ifelse(genome1 == genome2, genome1, "Intergenome")
        )]
        if(sourceSelection == "multiple"){ # to allow merging across multiple sources with different samples 
            cols <- names(x)
            x[, .SD, .SDcols = cols[!(cols %in% c(samples_$Sample_ID))]]           
        } else x
    }))
    stopSpinner(session)
    x
})
filteredJunctions <- reactive({ # junction subjected to page-level filters only
    x <- junctions() 
    startSpinner(session, message = "filtering junctions.")
    x <- svx_applyJunctionFilters(x, stepSettings)
    stopSpinner(session)
    x
})

#----------------------------------------------------------------------
# selected samples, stratified by genome
#----------------------------------------------------------------------
sampleGenomesTableData <- reactive({
    sourceIds <- sourceIds()
    req(sourceIds)
    samples <- samples()
    genomes <- genomes()
    do.call(rbind, lapply(sourceIds, function(sourceId_){
        genomes <- genomes[sourceId == sourceId_]
        isCompositeGenome <- genomes$isCompositeGenome[1]
        samples[
            sourceId == sourceId_, 
            .(
                Source_Genome = genomes$sourceGenome[1],
                Genome = if(isCompositeGenome) c(genomes$genome, "Intergenome") else genomes$genome, 
                isCompositeGenome = isCompositeGenome, 
                compositeDelimiter = genomes$compositeDelimiter[1],
                Sample = getSampleNames(sampleIds = Sample_ID)
            ), 
            by = .(sourceId, Project, Sample_ID)
        ]
    }))
})
sampleGenomesTable <- bufferedTableServer(
    "sampleGenomes",
    id,
    input,
    tableData = reactive( sampleGenomesTableData()[, .(Project, Sample, Source_Genome, Genome)] ),
    selection = sampleGenomesSelection,
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
)
sampleGenomes <- reactive({
    rows <- sampleGenomesTable$rows_selected()
    req(rows)
    sampleGenomesTableData()[rows]
})
sampleGenomeJunctions <- reactive({ # further apply the sampleGenome filters
    jxns <- filteredJunctions()
    sampleGenomes <- sampleGenomes()    
    startSpinner(session, message = "filtering junctions..")
    x <- jxns[sourceGenome %in% sampleGenomes$Source_Genome &
              Genome       %in% sampleGenomes$Genome]
    stopSpinner(session)
    x
})

#----------------------------------------------------------------------
# module return value
#----------------------------------------------------------------------
list(
    sourceIds = sourceIds,
    samples = samples,
    genomes = genomes,
    junctions = junctions,
    filteredJunctions = filteredJunctions,
    sampleGenomes = sampleGenomes,
    sampleGenomeJunctions = sampleGenomeJunctions
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
