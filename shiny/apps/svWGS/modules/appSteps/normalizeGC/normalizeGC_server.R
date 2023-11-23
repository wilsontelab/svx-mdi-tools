#----------------------------------------------------------------------
# server components for the normalizeGC appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
normalizeGCServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'normalizeGC'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
gridColors <- list(x = "#aaaaaa", y = "#aaaaaa")

#----------------------------------------------------------------------
# data package sources and source-level data objects derived from pipeline
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("source", selection = "single") 
sampleNames <- reactive({ # vector of the names of all co-analyzed samples
    sourceId <- sourceId()
    req(sourceId)
    source <- getSourceFromId(sourceId)
    req(source)
    source$manifest$Sample_ID
})
genomes <- reactive({ # data.table with genome,binSize,isComposite,compositeDelimiter
    sourceId <- sourceId()
    req(sourceId)
    readRDS(getSourceFilePath(sourceId, "genomesFile"))
})
chromosomes <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    readRDS(getSourceFilePath(sourceId, "chromosomesFile"))
})
bins <- reactive({ # contiguous genome bins of size --bin-size per genome, in one table, one column at end per co-analyzed sample
    sourceId <- sourceId()
    req(sourceId)    
    sampleNames <- sampleNames()
    req(sampleNames)
    startSpinner(session, message = "loading coverage")  
    binsFile <- getSourceFilePath(sourceId, "binsCoverageFile")   
    readRDS(binsFile)[, .SD, .SDcols = c("chrom","start","gc","excluded", sampleNames)] # not yet filtered by genome
})
junctions <- reactive({
    sourceId <- sourceId()
    req(sourceId)  
    x <- svx_loadJunctions(sourceId, svWGS_loadJunctions) 
    stopSpinner(session)
    x
})
nonExcludedBins <- reactive({
    bins()[excluded == 0]
})
getGenome <- function(genomeName) genomes()[genome == genomeName]
getGenomeBins <- function(genome, bins = NULL, includeExcluded = TRUE) {
    if(is.null(bins)) bins <- if(includeExcluded) bins() else nonExcludedBins()
    if(genome$isComposite) bins[endsWith(chrom, genome$genome)] else bins
}
getGenomeJxns <- function(genome, junctions = NULL) {
    if(is.null(junctions)) junctions <- junctions()
    if(genome$isComposite) junctions[
        endsWith(CHROM_1, genome$genome) & # allow translocations but suppress intergenomic
        endsWith(CHROM_2, genome$genome)
    ] else junctions
}

#----------------------------------------------------------------------
# source samples, stratified by composite genomes when applicable
#----------------------------------------------------------------------
sampleGenomeTableData <- reactive({
    as.data.table(expand.grid(
        Sample = sampleNames(), 
        Genome = genomes()$genome, 
        stringsAsFactors = FALSE
    ))[, ":="(
        key = paste(Sample, Genome),
        sourceId = sourceId()        
    )]
})
sampleGenomeTable <- bufferedTableServer(
    "sampleGenome",
    id,
    input,
    tableData = reactive( sampleGenomeTableData()[, .(Sample, Genome)] ),
    selection = 'single',
    options = list(
        paging = FALSE,
        searching = FALSE  
    )
)
sampleGenome <- reactive({
    row <- sampleGenomeTable$rows_selected()
    req(row)
    sampleGenomeTableData()[row]
})

#----------------------------------------------------------------------
# appStep outcomes, saved to disk since these are ~one-time analysis steps
# includes GC bias fit, chromosome-level data and junction fits, but not HMM
#----------------------------------------------------------------------
gcBiasFileName <- "gcBiasModels.rds"
gcBiasFile <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    expandSourceFilePath(sourceId, gcBiasFileName)
})
invalidateGcBiasModels <- reactiveVal(1)
getGcBiasModels <- function(gcBiasFile = NULL, sourceId = NULL){
    if(is.null(gcBiasFile)) gcBiasFile <- expandSourceFilePath(sourceId, gcBiasFileName)
    if(file.exists(gcBiasFile)) readRDS(gcBiasFile) else list()
}
gcBiasModels <- reactive({ # gc bias model from negative binomial are calculated synchronously
    invalidateGcBiasModels()
    getGcBiasModels( gcBiasFile = gcBiasFile() )
})
gcBiasModel <- reactive({
    gcBiasModels <- gcBiasModels()    
    sampleGenome <- sampleGenome()
    gcBiasModels[[sampleGenome$key]]
})
suppressGcOutliers <- function(nb, gc){
    allowedGC <- range(nb$model$fractionGC)    
    pmax(allowedGC[1], pmin(allowedGC[2], gc))
}

#----------------------------------------------------------------------
# interactive GC bias plot, selection cascades to solving negative binomial
#----------------------------------------------------------------------
gcPlotData <- function(){
    sampleGenome <- sampleGenome()
    bins <- getGenomeBins( getGenome(sampleGenome$Genome), includeExcluded = FALSE )
    startSpinner(session, message = paste("plotting", sampleGenome$key))
    d <- bins[, .SD, .SDcols = c("gc", sampleGenome$Sample)]   
    setnames(d, c("x", "y"))
    d[sample.int(.N, min(.N, 5000))]
}
gcOverplotData <- function(){
    gcBiasModel <- gcBiasModel()    
    if(!isTruthy(gcBiasModel)) return(NULL)    
    sampleGenome <- sampleGenome()
    startSpinner(session, message = paste("overplotting", sampleGenome$key))
    updateNumericInput(session, "ploidy", value = gcBiasModel$ploidy)
    updateNumericInput(session, "expectedSex", value = gcBiasModel$expectedSex)
    nb <- gcBiasModel$fit
    gc <- nb$model$fractionGC
    rpa <- predict(nb, gc, type = 'mu') # rpa = reads per allele
    data.table(
        x = c(gc, NA, gc, NA, gc, NA, gc), # draw predicted coverage traces for CN 1 to 4
        y = c(rpa * 1, NA, rpa * 2, NA, rpa * 3, NA, rpa * 4)
    )
}
gcBiasPlot <- interactiveScatterplotServer(
    "gcBiasPlot",
    plotData = reactive({ 
        x <- gcPlotData() 
        stopSpinner(session)
        x
    }),
    accelerate = TRUE,
    color = CONSTANTS$plotlyColors$blue,
    overplot = reactive({
        x <- gcOverplotData()
        stopSpinner(session)
        x
    }),
    overplotMode = "lines",
    overplotColor = CONSTANTS$plotlyColors$red,
    xtitle = "Fraction GC",
    xrange = c(0.2, 0.75),
    ytitle = "Read Depth",
    yrange = function(...) range_pos(..., foldIQR = 3),
    selectable = "lasso"
)    
observeEvent(gcBiasPlot$selected(), {
    showUserDialog(
        "Use this GC Bias Fit?", 
        tags$p(
            "If you are happy with your GC bias selection, click OK to run the next, slow actions."
        ), 
        tags$p(
            "If you are not happy, click Cancel and repeat the GC selection for the CN=ploidy group."
        ), 
        callback = function(parentInput) {
            removeModal()
            fitGCBiasFromSelected(gcBiasPlot$selected())
        },
        type = 'okCancel', 
        easyClose = FALSE, 
        fade = 250
    )
})

#----------------------------------------------------------------------
# fit negative binomial distribution to GC bias; normalize bins and SV junctions to the fit
#----------------------------------------------------------------------
fitGCBiasFromSelected <- function(selected){
    req(selected, nrow(selected) > 10)
    sampleGenome <- sampleGenome()
    genome <- getGenome(sampleGenome$Genome)    
    startSpinner(session, message = paste("fitting", sampleGenome$key))
    selected <- as.data.table(selected)
    gcBiasModels <- gcBiasModels()
    fit <- new_nbinomCountsGC(
        binCounts  = selected$y, 
        fractionGC = selected$x, 
        binCN = input$ploidy,
        method = 'cubic'
    )
    chromCN <- calculateChromCN(sampleGenome, genome, fit)
    jxnCN <- normalizeJxnCN(sampleGenome, genome, fit, chromCN)
    gcBiasModels[[sampleGenome$key]] <- list(
        sampleGenome    = sampleGenome,
        genome          = genome,
        ploidy          = input$ploidy,
        expectedSex     = input$expectedSex,
        fit             = fit,
        chromCN         = chromCN,
        jxnCN           = jxnCN
    )
    saveRDS(gcBiasModels, file = gcBiasFile())
    stopSpinner(session)
    invalidateGcBiasModels( invalidateGcBiasModels() + 1 ) 
    # target <- list(
    #     sampleGenome = sampleGenome,
    #     genome       = genome,
    #     ploidy       = input$ploidy,
    #     modelKey     = sampleGenome$key,
    #     model        = gcBiasModels[[sampleGenome$key]]$fit,
    #     binsFile     = getSourceFilePath(sampleGenome$sourceId, "binsCoverageFile"),
    #     svsFile      = getSourceFilePath(sampleGenome$sourceId, "structuralVariants"),
    #     jxnBinsFile  = getSourceFilePath(sampleGenome$sourceId, "junctionBinsCoverageFile"),
    #     jxnGcFile    = getSourceFilePath(sampleGenome$sourceId, "junctionGcFile"),
    #     entropy      = sample(1e8, 1)
    # ) 
}
calculateChromCN <- function(sampleGenome, genome, fit){
    bins <- getGenomeBins(genome, includeExcluded = FALSE)[
        !startsWith(toupper(chrom), "CHRM"), 
        .SD, 
        .SDcols = c("chrom","gc",sampleGenome$Sample)
    ]
    startSpinner(session, message = "profiling chromosomes")  
    rpa <- predict(fit, suppressGcOutliers(fit, bins$gc), type = 'mu')
    coverage <- bins[[sampleGenome$Sample]]
    x <- bins[, 
        cn := coverage / rpa
    ][, 
        .(predominantCN = round(peakValue(cn), 0)), 
        by = .(chrom)
    ]
    x[, expectedCN := input$ploidy]  
    if(any(startsWith(toupper(x$chrom), "CHRY"))){ # do nothing if not an XY sex chromosome system
        expectedSex <- if(input$expectedSex == "auto"){ # in auto mode, determine expecteSex from chrY
            chrYCN <- x[startsWith(toupper(chrom), "CHRY"), predominantCN]
            if(chrYCN == 0) "XX" else "XY"
        } else input$expectedSex
        if(expectedSex == "XX"){
            x[startsWith(toupper(chrom), "CHRX"), expectedCN := input$ploidy]  
            x[startsWith(toupper(chrom), "CHRY"), expectedCN := 0]  
        } else if(expectedSex == "XY"){
            x[startsWith(toupper(chrom), "CHRX"), expectedCN := input$ploidy - 1]  
            x[startsWith(toupper(chrom), "CHRY"), expectedCN := 1]  
        }
    } 
    setkey(x, chrom)
    x
}
getNormalizedSpanCn <- function(nb, gc, coverage, overlap = NULL){
    rpa <- predict(nb, suppressGcOutliers(nb, gc), type = 'mu')
    cn <- coverage / rpa
    if(is.null(overlap)) mean(cn) else weighted.mean(cn, overlap)
}
normalizeJxnCN <- function(sampleGenome, genome, fit, chromCN){
    startSpinner(session, message = "normalizing junctions") 
    nonSingletonJxns <- getGenomeJxns(genome)[N_TOTAL - N_OUTER_CLIPS > 1] # never analyze singleton junctions, they aren't expected to impact CN
    nonSingletonJxns[, svId := SV_ID]
    setkey(nonSingletonJxns, svId)
    sampleBins <- getGenomeBins(genome, includeExcluded = FALSE)
    jxnBinsFile <- getSourceFilePath(sampleGenome$sourceId, "junctionBinsCoverageFile")
    jxnGcFile   <- getSourceFilePath(sampleGenome$sourceId, "junctionGcFile")
    sampleJxnBins <- readRDS(jxnBinsFile)[sample == sampleGenome$Sample & svId  %in% nonSingletonJxns[, svId]]
    sampleJxnGc   <- readRDS(jxnGcFile)[  sample == sampleGenome$Sample & SV_ID %in% nonSingletonJxns[, svId]][, svId := SV_ID]
    setkey(sampleJxnBins, svId) # svId for data merging
    setkey(sampleBins, chrom)
    setindex(sampleBins, start)
    setkey(sampleJxnGc, svId)
    x <- sampleJxnBins[, {
        gc       <- as.numeric(strsplit(gc, ",")[[1]])
        coverage <- as.numeric(strsplit(coverage, ",")[[1]])
        overlap  <- as.integer(strsplit(overlap, ",")[[1]])
        .( 
            cn  = getNormalizedSpanCn(fit, gc, coverage, overlap)
        )
    }, by = .(svId, spanType)] %>%
    dcast(svId ~ spanType, fun.aggregate = sum, value.var = "cn") %>%
    merge(
        nonSingletonJxns[, .(
            svId,
            CHROM_1,
            CHROM_2,
            edgeType,
            SV_SIZE,
            N_SAMPLES,
            N_TOTAL,
            N_SPLITS
        )],
        by = "svId",
        all.y = TRUE
    ) %>%
    merge(
        sampleJxnGc,
        by = "svId",
        all.x = TRUE
    )
    setkey(x, svId)    
    x[, junctionCN := getNormalizedSpanCn(fit, gc, N_TOTAL), by = .(svId)]   
    x[, .(
        svId,
        CHROM_1,
        CHROM_2,
        edgeType,
        SV_SIZE,
        N_SAMPLES,
        N_TOTAL,
        N_SPLITS,
        gc,
        junctionCN,
        outerFlankCN = outerFlank,
        innerFlankCN = innerFlank,
        eventSpanCN  = eventSpan,
        flankCNC     = innerFlank - outerFlank,
        eventCNC     = eventSpan  - outerFlank,
        chromCNC     = eventSpan - chromCN[CHROM_1, predominantCN]
    )] 
}

#----------------------------------------------------------------------
# interactive chromosome bin density plot
#----------------------------------------------------------------------
observeEvent(chromosomes(), {
    chromosomes <- chromosomes()
    req(chromosomes)
    updateRadioButtons(session, "densityChrom", 
                       choices = c("genome", unique(chromosomes$canonicalChroms)), 
                       selected = "genome", inline = FALSE)
})
normalizedBins <- reactive({
    gcBiasModel <- gcBiasModel()  
    req(gcBiasModel)      
    sampleGenome <- sampleGenome()
    startSpinner(session, message = paste("calculating density", sampleGenome$key))
    nb <- gcBiasModel$fit
    bins <- getGenomeBins(gcBiasModel$genome)
    rpa <- predict(nb, suppressGcOutliers(nb, bins$gc), type = 'mu')
    x <- bins[, cn := bins[[sampleGenome$Sample]] / rpa] 
    stopSpinner(session)
    x
})
chromDensityPlotData <- reactive({
    bins <- normalizedBins()
    if(input$densityChrom != "genome") bins <- bins[chrom == input$densityChrom]
    as.data.table(density(bins$cn, na.rm = TRUE)[c("x","y")])
})
chromDensityPlotVLines <- reactive({
    if(input$densityChrom == "genome") return(input$ploidy)
    gcBiasModels <- gcBiasModels()    
    sampleGenome <- sampleGenome()
    gcBiasModel <- gcBiasModels[[sampleGenome$key]]
    req(gcBiasModel)
    gcBiasModel$chromCN[chrom == input$densityChrom, expectedCN]
})
chromDensityPlot <- interactiveScatterplotServer(
    "chromDensityPlot",
    plotData = chromDensityPlotData,
    mode = "lines",
    accelerate = TRUE,
    color = CONSTANTS$plotlyColors$blue,
    xtitle = "Copy Number",
    xrange = c(0, 4.5),
    ytitle = "Density",
    vLines = chromDensityPlotVLines,
    grid = gridColors
)  

# #----------------------------------------------------------------------
# # hmmCnvsFileName <- "hmmCnvs.rds"
# # hmmCnvsFile <- reactive({
# #     sourceId <- sourceId()
# #     req(sourceId)
# #     expandSourceFilePath(sourceId, hmmCnvsFileName)
# # })
# # invalidateHmmCnvs <- reactiveVal(1)
# # getHmmCnvs <- function(){
# #     hmmCnvsFile <- hmmCnvsFile()
# #     if(file.exists(hmmCnvsFile)) readRDS(hmmCnvsFile) else list()    
# # }
# # updateHmmCnvs <- function(sample, data){
# #     hmmCnvs <- getHmmCnvs()
# #     hmmCnvs[[sample]] <- data
# #     saveRDS(hmmCnvs, hmmCnvsFile())
# #     invalidateHmmCnvs( invalidateHmmCnvs() + 1 )
# # }
# # hmmCnvs <- reactive({ # bin normalized CN and HMM CNVs calls are calculated asynchronoulsy and stored separately
# #     invalidateHmmCnvs()
# #     getHmmCnvs()
# # })
# #----------------------------------------------------------------------
# # cascading, asynchrononous normalization functions
# #------------------------------------------------------------ ----------
# # solveCnvHMM <- function(target, chromCN){
# #     # bins <- readRDS(target$collapsedBinsFile)[, .SD, .SDcols = c("chrom","start","gc",target$sample)]
# #     bins <- readRDS(target$binsFile)[, .SD, .SDcols = c("chrom","start","gc",target$sample)]

# #     setnames(bins, c("chrom","start","gc","coverage"))
# #     binSize <- bins[1:2, diff(start)]
    
# #     # do.call(rbind, lapply(unique(bins$chrom), function(chrom_){
# #     do.call(rbind, lapply("chr12", function(chrom_){

# #         bins <- bins[chrom == chrom_]
# #         cn <- viterbi(
# #             target$model, 
# #             binCounts = bins$coverage, 
# #             fractionGC = suppressGcOutliers(target$model, bins$gc), 
# #             maxCN = 5,
# #             asRle = TRUE
# #         )$cn
# #         nSegments <- length(cn$lengths)
# #         if(nSegments == 1 && cn$values[1] == 2) return(NULL)
# #         endBins <- cumsum(cn$lengths)
# #         startBins <- c(1, endBins + 1)[1:nSegments]
# #         d <- do.call(rbind, lapply(1:nSegments, function(i){
# #             if(cn$values[i] == 2) return(NULL)
# #             bins[startBins[i]:endBins[i], .(
# #                 chrom = chrom_,
# #                 start = as.integer(start[1]), # BED format
# #                 end   = as.integer(start[.N] + binSize),
# #                 gc    = mean(gc),
# #                 count = sum(coverage),
# #                 hmmCn = cn$values[i],
# #                 cn    = getNormalizedSpanCn(target$model, gc, coverage)
# #             )]
# #         }))
# #         d[, ":="(
# #             svId       = paste("cnv", 1:.N, sep = ":"),
# #             SV_SIZE    = end - start,
# #             N_SAMPLES  = 1,
# #             nInstances = NA_integer_,
# #             cnc        = cn - chromCN[d$chrom, expectedCN] # event copy number relative to Mendelian expectations
# #         )]
# #         d[, ":="(
# #             edgeType   = ifelse(cnc > 0, "U", "D")
# #         )]
# #     }))
# # }
# # cnvHmmResults <- reactiveVal()
# # triggerCnvHmm <- function(target){
# #     req(target, target$sample)
# #     updateHmmCnvs(target$sample, NULL) # clear existing data and plot
# #     mdi_async(
# #         name = paste("hmm", target$sample),
# #         header = TRUE, 
# #         autoClear = 2000,
# #         taskFn = function(target) {
# #             chromCN <- calculateChromCN(target)
# #             list(
# #                 sample = target$sample,
# #                 chromCN = chromCN, # yes, repeat this here, since we don't know the order the user will do the two GC selections
# #                 cnvHmm = solveCnvHMM(target, chromCN)
# #             )
# #         },
# #         target = target,
# #         reactiveVal = cnvHmmResults
# #     )
# # }
# # observeEvent(cnvHmmResults(), {
# #     results <- cnvHmmResults()
# #     req(results, results$value, results$success)
# #     updateSampleChromCN(results$value$sample, results$value$chromCN)
# #     updateHmmCnvs(results$value$sample, results$value$cnvHmm)
# # })

#----------------------------------------------------------------------
# external utilities for browser support, etc.
#----------------------------------------------------------------------
getGcBiasModels_externalCall <- function(sourceId){
    gcSourceId <- tryCatch(sourceId(), error = function(e) NULL)
    if(isTruthy(gcSourceId) && gcSourceId == sourceId) gcBiasModels()
    else getGcBiasModels(sourceId = sourceId)
}
getBinNormalizedCN <- function(sourceId, sampleName, reference, coord, gc, coverage){ # for plotting GC-normalized copy number at bin level
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    if(isCompositeGenome(reference)){
        modelKey <- paste(sampleName, coord$chromosome) # in case one of multiple genomes is selected as chrom
        gcBiasModel <- gcBiasModels[[modelKey]]
        if(!isTruthy(gcBiasModel)){
            delim <- getCustomCompositeDelimiter(reference$metadata)
            if(grepl(delim, coord$chromosome)){
                modelKey <- paste(sampleName, strsplit(coord$chromosome, delim)[[1]][2])
                gcBiasModel <- gcBiasModels[[modelKey]]
            }
        }
    } else {
        modelKey <- paste(sampleName, reference$genome$genome)
        gcBiasModel <- gcBiasModels[[modelKey]]
    }
    if(!isTruthy(gcBiasModel)) return(NULL)
    rpa <- predict(gcBiasModel$fit, suppressGcOutliers(gcBiasModel$fit, gc), type = 'mu')
    coverage / rpa
}
getCnvJxnsNormalizedCN <- function(sourceId){
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    sessionCache$get(
        'getCnvJxnsNormalizedCN', 
        keyObject = list(sourceId = sourceId, gcBiasModels = gcBiasModels), # thus, this cached object updates as gc bias is refit
        permanent = FALSE, 
        from = "ram",
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "setting junction CN")
            x <- do.call(rbind, lapply(names(gcBiasModels), function(x) {
                gcBiasModels[[x]]$jxnCN[, .(svId, outerFlankCN, innerFlankCN, flankCNC, junctionCN)]
            }))
            x[, .(
                outerFlankCN = mean(outerFlankCN, na.rm = TRUE), 
                innerFlankCN = mean(innerFlankCN, na.rm = TRUE), 
                flankCNC     = mean(flankCNC, na.rm = TRUE),
                junctionCN   = mean(junctionCN, na.rm = TRUE)
            ), by = .(svId)][, SV_ID := svId]            
        }
    )
}

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    # updateTextInput(session, 'xxx', value = bm$outcomes$xxx)
    # xxx <- bm$outcomes$xxx
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = list(
    ),
    # isReady = reactive({ getStepReadiness(options$source, gcBiasModels()) }),
    getGcBiasModels = getGcBiasModels_externalCall,
    getBinNormalizedCN = getBinNormalizedCN,
    getCnvJxnsNormalizedCN = getCnvJxnsNormalizedCN,
    # getNormalizedHmmCnvs = getNormalizedHmmCnvs,
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
