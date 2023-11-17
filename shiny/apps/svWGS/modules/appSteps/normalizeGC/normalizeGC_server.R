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

#----------------------------------------------------------------------
# sample sources and data derived from pipeline
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
bins <- reactive({ # contiguous genome bins of size --bin-size per genome, in one table, one column at end per co-analyzed sample
    sourceId <- sourceId()
    req(sourceId)    
    sampleNames <- sampleNames()
    req(sampleNames)
    startSpinner(session, message = "loading coverage")  
    binsFile <- getSourceFilePath(sourceId, "coverageFile")   
    readRDS(binsFile)[, .SD, .SDcols = c("chrom","start","gc", sampleNames)] # not yet filtered by genome
})
getGenomeBins <- function(genome, bins = NULL) {
    if(is.null(bins)) bins <- bins()
    if(genome$isComposite) bins[endsWith(chrom, genome$genome)] else bins
}
getCnvJxnBins <- function(genome, junctions) {
    if(genome$isComposite) junctions[endsWith(CHROM_1, genome$genome)] else junctions
}

#----------------------------------------------------------------------
# appStep outcomes, saved to disk since these are ~one-time analysis steps
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
suppressGcOutliers <- function(nb, gc){
    allowedGC <- range(nb$model$fractionGC)    
    pmax(allowedGC[1], pmin(allowedGC[2], gc))
}
#----------------------------------------------------------------------
chromCnFileName <- "chromCn.rds"
getChromCnFile <- function(sourceId, genome){
    req(sourceId, genome)
    expandSourceFilePath(sourceId, paste(genome$genome, chromCnFileName, sep = "."))
}
getChromCN <- function(sourceId, genome){
    chromCnFile <- getChromCnFile(sourceId, genome)
    if(file.exists(chromCnFile)) readRDS(chromCnFile) else list()
}
updateSampleChromCN <- function(sample, genome, data){
    sourceId <- sourceId()
    chromCN <- getChromCN(sourceId, genome)
    dataKey <- paste(sample, genome$genome)
    chromCN[[dataKey]] <- data
    saveRDS(chromCN, getChromCnFile(sourceId, genome))
}
#----------------------------------------------------------------------
normalizedCnvsFileName <- "normalizedCnvs.rds"
getNormalizedCnvsFile <- function(sourceId, genome){
    req(sourceId, genome)
    expandSourceFilePath(sourceId, paste(genome$genome, normalizedCnvsFileName, sep = "."))
}
invalidateNormalizedCnvs <- reactiveVal(1)
getNormalizedCnvs <- function(sourceId, genome){
    normalizedCnvsFile <- getNormalizedCnvsFile(sourceId, genome)
    if(file.exists(normalizedCnvsFile)) readRDS(normalizedCnvsFile) else list()    
}
updateSampleNormalizedCnvs <- function(sample, genome, data){
    sourceId <- sourceId()
    normalizedCnvs <- getNormalizedCnvs(sourceId, genome)
    dataKey <- paste(sample, genome$genome)
    normalizedCnvs[[dataKey]] <- data
    saveRDS(normalizedCnvs, getNormalizedCnvsFile(sourceId, genome))
    invalidateNormalizedCnvs( invalidateNormalizedCnvs() + 1 )
}
#----------------------------------------------------------------------
hmmCnvsFileName <- "hmmCnvs.rds"
hmmCnvsFile <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    expandSourceFilePath(sourceId, hmmCnvsFileName)
})
invalidateHmmCnvs <- reactiveVal(1)
getHmmCnvs <- function(){
    hmmCnvsFile <- hmmCnvsFile()
    if(file.exists(hmmCnvsFile)) readRDS(hmmCnvsFile) else list()    
}
updateHmmCnvs <- function(sample, data){
    hmmCnvs <- getHmmCnvs()
    hmmCnvs[[sample]] <- data
    saveRDS(hmmCnvs, hmmCnvsFile())
    invalidateHmmCnvs( invalidateHmmCnvs() + 1 )
}
hmmCnvs <- reactive({ # bin normalized CN and HMM CNVs calls are calculated asynchronoulsy and stored separately
    invalidateHmmCnvs()
    getHmmCnvs()
})

#----------------------------------------------------------------------
# interactive GC bias plot, selection cascades to solving negative binomial
#----------------------------------------------------------------------
gcPlotObservers <- list()
gcPlotData <- function(sampleName, genome){
    bins <- getGenomeBins(genome)
    startSpinner(session, message = paste("plotting", sampleName, genome$genome))
    d <- bins[, .SD, .SDcols = c("gc", sampleName)]   
    setnames(d, c("x", "y"))
    d[sample.int(.N, min(.N, 5000))]
}
gcOverplotData <- function(sourceId, sampleName, genome){
    gcBiasModels <- gcBiasModels()
    modelKey <- paste(sampleName, genome$genome)
    fit <- gcBiasModels[[modelKey]]
    if(!isTruthy(fit)) return(NULL)
    startSpinner(session, message = paste("overplotting", sampleName, genome$genome))
    gc <- fit$model$fractionGC
    rpa <- predict(fit, gc, type = 'mu') # rpa = reads per allele
    data.table(
        x = c(gc, NA, gc, NA, gc), # draw predicted coverage traces for CN 1 to 3
        y = c(rpa * 1, NA, rpa * 2, NA, rpa * 3)
    )
}
gcPlotServer <- function(id, sourceId, sampleName, genome){
    modelKey <- paste(sampleName, genome$genome)
    plot <- interactiveScatterplotServer(
        id,
        plotData = reactive({ 
            x <- gcPlotData(sampleName, genome) 
            stopSpinner(session)
            x
        }),
        accelerate = TRUE,
        color = CONSTANTS$plotlyColors$blue,
        overplot = reactive({
            x <- gcOverplotData(sourceId, sampleName, genome)
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
    if(!is.null(gcPlotObservers[[modelKey]])) gcPlotObservers[[modelKey]]$destroy()
    gcPlotObservers[[modelKey]] <<- observeEvent(plot$selected(), {
        showUserDialog(
            "Use this GC Bias Fit?", 
            tags$p(
                "If you are happy with your GC bias selection, click OK to run the next, slow actions."
            ), 
            tags$p(
                "If you are not happy, click Cancel and repeat the GC selection for the CN=ploidy group."
            ), 
            callback = function(parentInput) fitGCBiasFromSelected(sourceId, sampleName, genome, plot$selected()),
            type = 'okCancel', 
            easyClose = FALSE, 
            fade = 250
        )
    })
}

#----------------------------------------------------------------------
# interactive junction plot, correlating predicted CNV sizes to normalized CN
#----------------------------------------------------------------------
jxnPlotData <- function(sourceId, sampleName, genome, jxnType){
    normalizedCnvs <- getNormalizedCnvs(sourceId, genome)
    dataKey <- paste(sampleName, genome$genome)
    req(normalizedCnvs, normalizedCnvs[[dataKey]])
    d <- normalizedCnvs[[dataKey]][N_SAMPLES == 1 & JXN_TYPE == jxnType, .( # therefore, CNV junctions unique to this sample
        x = log10(SV_SIZE),
        y = flankCNC
    )]
}
jxnPlotServer <- function(id, sourceId, sampleName, genome){
    plot <- interactiveScatterplotServer(
        id,
        plotData = reactive({ 
            x <- jxnPlotData(sourceId, sampleName, genome, "L") 
            stopSpinner(session)
            x
        }),
        accelerate = TRUE,
        pointSize = 5,
        overplot = reactive({
            jxnPlotData(sourceId, sampleName, genome, "D")
        }),
        overplotPointSize = 5,
        xtitle = "log10(CNV Size)",
        xrange = c(2.75, 7),
        ytitle = "Copy Number Change",
        yrange = c(-2.5, 2.5),
        grid  = list(x ="grey", y = "grey")
    ) 
}

#----------------------------------------------------------------------
# interactive HMM CNV plot, correlating observed CNV sizes to normalized CN
#----------------------------------------------------------------------
cnvPlotData <- function(sampleName, jxnType){
    hmmCnvs <- hmmCnvs()
    req(hmmCnvs, hmmCnvs[[sampleName]])
    d <- hmmCnvs[[sampleName]][JXN_TYPE == jxnType, .( # therefore, CNV junctions unique to this sample
        x = log10(SV_SIZE),
        y = cnc
    )]
}
cnvPlotServer <- function(id, sourceId, sampleName, genome){

return(NULL)

    plot <- interactiveScatterplotServer(
        id,
        plotData = reactive({ 
            x <- cnvPlotData(sampleName, "L") 
            stopSpinner(session)
            x
        }),
        accelerate = TRUE,
        pointSize = 5,
        overplot = reactive({
            cnvPlotData(sampleName, "D")
        }),
        overplotPointSize = 5,
        xtitle = "log10(CNV Size)",
        xrange = c(4.75, log10(250e6)),
        ytitle = "Copy Number Change",
        yrange = c(-2.5, 2.5),
        grid  = list(x ="grey", y = "grey")
    ) 
}

#----------------------------------------------------------------------
# one normalization box per sample, with GC and CNV plots
#----------------------------------------------------------------------
normalizePlotUI <- function(id, title){
    box(
        width = 4,
        title = title,
        interactiveScatterplotUI(session$ns(id), height = '400px')
    )        
}
output$samples <- renderUI({
    sourceId <- sourceId()
    sampleNames <- sampleNames()
    genomes <- genomes()
    req(sourceId, sampleNames, genomes)
    lapply(sampleNames, function(sampleName){
        lapply(1:nrow(genomes), function(genomeI){
            genome <- genomes[genomeI]
            gcBiasId  <- paste("gcBias",  sampleName, genome, sep = "_")
            jxnPlotId <- paste("jxnPlot", sampleName, genome, sep = "_")
            cnvPlotId <- paste("cnvPlot", sampleName, genome, sep = "_")
            gcPlotServer( gcBiasId,  sourceId, sampleName, genome)
            jxnPlotServer(jxnPlotId, sourceId, sampleName, genome)
            cnvPlotServer(cnvPlotId, sourceId, sampleName, genome)
            fluidRow(box(
                width = 12,
                title = paste0(sampleName, " (", genome$genome, ")"),
                solidHeader = TRUE,
                status = "primary",
                normalizePlotUI(gcBiasId,  paste0("GC Bias Plot (", as.integer(genome$binSize / 1e3), "kb, downsampled)")),
                normalizePlotUI(jxnPlotId, "Unique Del/Dup Junctions"),
                normalizePlotUI(cnvPlotId, "All HMM CNVs")
            ))
        })
    })
})

#----------------------------------------------------------------------
# cascading, asynchrononous normalization functions
#----------------------------------------------------------------------
fitGCBiasFromSelected <- function(sourceId, sampleName, genome, selected, triggerFn){
    req(selected, nrow(selected) > 10)
    startSpinner(session, message = paste("fitting", sampleName, genome$genome))
    selected <- as.data.table(selected)
    gcBiasModels <- gcBiasModels()
    modelKey <- paste(sampleName, genome$genome)
    gcBiasModels[[modelKey]] <- new_nbinomCountsGC(binCounts = selected$y, fractionGC = selected$x, method = 'cubic')
    saveRDS(gcBiasModels, file = gcBiasFile())
    stopSpinner(session)
    invalidateGcBiasModels( invalidateGcBiasModels() + 1 ) 
    target <- list(
        sourceId = sourceId,
        sample   = sampleName, 
        genome   = genome,
        modelKey = modelKey,
        model    = gcBiasModels[[modelKey]],
        binsFile    = getSourceFilePath(sourceId, "coverageFile"),
        svsFile     = getSourceFilePath(sourceId, "structuralVariants"),
        cnvBinsFile = getSourceFilePath(sourceId, "cnvCoverageFile"),
        entropy  = sample(1e8, 1)
    ) 
    triggerCnvNormalization(target)
}
#----------------------------------------------------------------------
calculateChromCN <- function(target){
    bins <- getGenomeBins(target$genome, readRDS(target$binsFile))[
        !startsWith(toupper(chrom), "CHRM"), 
        .SD, 
        .SDcols = c("chrom","gc",target$sample)
    ]
    rpa <- predict(target$model, suppressGcOutliers(target$model, bins$gc), type = 'mu')
    coverage <- bins[[target$sample]]
    hasChrY <- any(startsWith(toupper(bins$chrom), "CHRY"))
    x <- bins[, cn := coverage / rpa][, 
        {
            predominantCN <- as.integer(peakValue(cn))
            if(hasChrY && any(startsWith(toupper(chrom), c("CHRX","CHRY")))) .(
                expectedCN    = predominantCN,
                predominantCN = predominantCN
            ) else .(
                expectedCN = 2L, # TODO: need ploidy input/outcome (not in pipeline, as data may change that answer)
                predominantCN = predominantCN
            )
        },
        by = .(chrom)
    ]
    setkey(x, chrom)
    x
}
getNormalizedSpanCn <- function(nb, gc, coverage, overlap = NULL){
    rpa <- predict(nb, suppressGcOutliers(nb, gc), type = 'mu')
    cn <- coverage / rpa
    if(is.null(overlap)) mean(cn) else weighted.mean(cn, overlap)
}
normalizeCnvCN <- function(target, chromCN){
    nonSingletonCnvJxns <- getCnvJxnBins(target$genome, readRDS(target$svsFile))[ # never analyze singleton junctions, they aren't expected to impact CN
        JXN_TYPE %in% c("D", "L"), # D == duplication, L = deletion
        .(
            svId = SV_ID,
            CHROM_1,
            minPos = pmin(POS_1, POS_2),
            maxPos = pmax(POS_1, POS_2),
            JXN_TYPE,
            SV_SIZE,
            N_SAMPLES,
            nInstances = N_TOTAL - N_OUTER_CLIPS
        )
    ][nInstances > 1]
    setkey(nonSingletonCnvJxns, svId)
    sampleBins    <- getGenomeBins(target$genome, readRDS(target$binsFile))
    sampleCnvBins <- readRDS(target$cnvBinsFile)[sample == target$sample & svId %in% nonSingletonCnvJxns[, svId]]
    setkey(sampleCnvBins, svId) # svId matches cnvJxns for data merging
    setkey(sampleBins, chrom)
    setindex(sampleBins, start)
    setindex(sampleBins, end)
    x <- merge(
        sampleCnvBins[, {
            gc       <- as.numeric(strsplit(gc, ",")[[1]])
            coverage <- as.numeric(strsplit(coverage, ",")[[1]])
            overlap  <- as.integer(strsplit(overlap, ",")[[1]])
            .( 
                cn  = getNormalizedSpanCn(target$model, gc, coverage, overlap)
            )
        }, by = .(svId)],
        nonSingletonCnvJxns,
        by = "svId",
        all.x = TRUE
    )
    binSize <- target$genome$binSize 
    elevenBinsSize <- binSize * 11
    x[, flankCN := sampleBins[CHROM_1][
            end   %between% c(minPos - elevenBinsSize, minPos - binSize) | 
            start %between% c(maxPos + binSize, maxPos + elevenBinsSize), # thus nominally 20 cleanly flanking bins
            getNormalizedSpanCn(target$model, gc, .SD[[target$sample]])
    ], by = .(svId)]
    x[, ":="(
        flankCNC = cn - flankCN,                        # event copy number relative to the DNA just outside of the junction span
        chromCNC = cn - chromCN[CHROM_1, predominantCN] # event copy number relative to the chromosome at large
    )]
}
cnvNormalizationResults <- reactiveVal()
triggerCnvNormalization <- function(target){
    req(target, target$sample)
    updateSampleNormalizedCnvs(target$sample, target$genome, NULL) # clear existing data and plot
    mdi_async(
        name = paste("normalize", target$sample, target$genome$genome),
        header = TRUE, 
        autoClear = 2000,
        taskFn = function(target) {
            chromCN <- calculateChromCN(target)
            list(
                target = target,
                chromCN = chromCN,
                cnvCN   = normalizeCnvCN(target, chromCN)
            )
        },
        target = target,
        reactiveVal = cnvNormalizationResults
    )
}
observeEvent(cnvNormalizationResults(), {
    results <- cnvNormalizationResults()
    req(results, results$value, results$success)
    x <- results$value
    updateSampleChromCN(x$target$sample, x$target$genome, x$chromCN)
    updateSampleNormalizedCnvs(x$target$sample, x$target$genome, x$cnvCN)
})
#----------------------------------------------------------------------
solveCnvHMM <- function(target, chromCN){
    # bins <- readRDS(target$collapsedBinsFile)[, .SD, .SDcols = c("chrom","start","gc",target$sample)]
    bins <- readRDS(target$binsFile)[, .SD, .SDcols = c("chrom","start","gc",target$sample)]

    setnames(bins, c("chrom","start","gc","coverage"))
    binSize <- bins[1:2, diff(start)]
    
    # do.call(rbind, lapply(unique(bins$chrom), function(chrom_){
    do.call(rbind, lapply("chr12", function(chrom_){

        bins <- bins[chrom == chrom_]
        cn <- viterbi(
            target$model, 
            binCounts = bins$coverage, 
            fractionGC = suppressGcOutliers(target$model, bins$gc), 
            maxCN = 5,
            asRle = TRUE
        )$cn
        nSegments <- length(cn$lengths)
        if(nSegments == 1 && cn$values[1] == 2) return(NULL)
        endBins <- cumsum(cn$lengths)
        startBins <- c(1, endBins + 1)[1:nSegments]
        d <- do.call(rbind, lapply(1:nSegments, function(i){
            if(cn$values[i] == 2) return(NULL)
            bins[startBins[i]:endBins[i], .(
                chrom = chrom_,
                start = as.integer(start[1]), # BED format
                end   = as.integer(start[.N] + binSize),
                gc    = mean(gc),
                count = sum(coverage),
                hmmCn = cn$values[i],
                cn    = getNormalizedSpanCn(target$model, gc, coverage)
            )]
        }))
        d[, ":="(
            svId       = paste("cnv", 1:.N, sep = ":"),
            SV_SIZE    = end - start,
            N_SAMPLES  = 1,
            nInstances = NA_integer_,
            cnc        = cn - chromCN[d$chrom, expectedCN] # event copy number relative to Mendelian expectations
        )]
        d[, ":="(
            JXN_TYPE   = ifelse(cnc > 0, "D", "L")
        )]
    }))
}
cnvHmmResults <- reactiveVal()
triggerCnvHmm <- function(target){
    req(target, target$sample)
    updateHmmCnvs(target$sample, NULL) # clear existing data and plot
    mdi_async(
        name = paste("hmm", target$sample),
        header = TRUE, 
        autoClear = 2000,
        taskFn = function(target) {
            chromCN <- calculateChromCN(target)
            list(
                sample = target$sample,
                chromCN = chromCN, # yes, repeat this here, since we don't know the order the user will do the two GC selections
                cnvHmm = solveCnvHMM(target, chromCN)
            )
        },
        target = target,
        reactiveVal = cnvHmmResults
    )
}
observeEvent(cnvHmmResults(), {
    results <- cnvHmmResults()
    req(results, results$value, results$success)
    updateSampleChromCN(results$value$sample, results$value$chromCN)
    updateHmmCnvs(results$value$sample, results$value$cnvHmm)
})

#----------------------------------------------------------------------
# external utilities for browser support, etc.
#----------------------------------------------------------------------
getGcBiasModels_externalCall <- function(sourceId){
    gcSourceId <- tryCatch(sourceId(), error = function(e) NULL)
    if(isTruthy(gcSourceId) && gcSourceId == sourceId) gcBiasModels()
    else getGcBiasModels(sourceId = sourceId)
}
getBinNormalizedCN <- function(sourceId, sampleName, genomeName, gc, coverage){ # for plotting GC-normalized copy number at bin level
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    modelKey <- paste(sampleName, genomeName)
    fit <- gcBiasModels[[modelKey]]
    if(!isTruthy(fit)) return(NULL)
    rpa <- predict(fit, suppressGcOutliers(fit, gc), type = 'mu')
    coverage / rpa
}
getCnvJxnsNormalizedCN <- function(sourceId, genomeName){
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    sessionCache$get(
        'getCnvJxnsNormalizedCN', 
        keyObject = list(sourceId = sourceId, genomeName = genomeName, gcBiasModels = gcBiasModels), # thus, this cached object updates as gc bias is refit
        permanent = FALSE, 
        from = "ram",
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "setting junction CN")
            normalizedCnvsFile <- expandSourceFilePath(sourceId, paste(genomeName, normalizedCnvsFileName, sep = "."))
            if(!file.exists(normalizedCnvsFile)) return(NA)
            normalizedCnvs <- readRDS(normalizedCnvsFile)
            samples <- names(normalizedCnvs)
            dt <- Reduce(
                function(dt1, dt2) merge(dt1, dt2, by = "svId", all.x = TRUE, all.y = TRUE),
                lapply(samples, function(sample) normalizedCnvs[[sample]][, .(svId, flankCNC)])
            )
            names(dt) <- c("SV_ID", samples)
            m <- dt[, .SD, .SDcols = samples]
            list(
                dt = dt[, ":="(
                    meanCNC = apply(m, 1, mean, na.rm = TRUE),
                    maxCNC  = apply(m, 1, function(v) max(abs(v), 0, na.rm = TRUE))
                )],
                fittedSamples = samples
            )
        }
    )
}
getNormalizedHmmCnvs <- function(sourceId){
    gcBiasModels <- getGcBiasModels_externalCall(sourceId)
    sessionCache$get(
        'getNormalizedHmmCnvs', 
        keyObject = list(sourceId = sourceId, gcBiasModels = gcBiasModels), # thus, this cached object updates as gc bias is refit
        permanent = FALSE, 
        from = "ram",
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "extracting HMM CNVs")
            hmmCnvsFile <- expandSourceFilePath(sourceId, hmmCnvsFileName)
            if(!file.exists(hmmCnvsFile)) return(NA)
            readRDS(hmmCnvsFile)

            # hmmCnvs <- readRDS(hmmCnvsFile)
            # samples <- names(hmmCnvs)


            # dt <- Reduce(
            #     function(dt1, dt2) merge(dt1, dt2, by = "svId", all.x = TRUE, all.y = TRUE),
            #     lapply(samples, function(sample) normalizedCnvs[[sample]][, .(svId, flankCNC)])
            # )
            # names(dt) <- c("SV_ID", samples)
            # m <- dt[, .SD, .SDcols = samples]
            # list(
            #     dt = dt[, ":="(
            #         meanCNC = apply(m, 1, mean, na.rm = TRUE),
            #         maxCNC  = apply(m, 1, function(v) max(abs(v), na.rm = TRUE))
            #     )],
            #     fittedSamples = samples
            # )
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
    outcomes = list(),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    getBinNormalizedCN = getBinNormalizedCN,
    getCnvJxnsNormalizedCN = getCnvJxnsNormalizedCN,
    getNormalizedHmmCnvs = getNormalizedHmmCnvs,
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
