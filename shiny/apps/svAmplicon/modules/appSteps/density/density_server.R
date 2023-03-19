#----------------------------------------------------------------------
# server components for the density appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
densityServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'density'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)
CONSTANTS$svTypes <- list(
    Translocation = "T", # edge/junction types (might be several per source molecule)
    Inversion     = "V",
    Duplication   = "U",
    Deletion      = "D",
    Unknown       = "?",
    Proper        = "P",
    Insertion     = "I", 
    MergeFailure  = "M"
)

#----------------------------------------------------------------------
# get the samples and data to plot
#----------------------------------------------------------------------
sampleSet <- sampleSetServer(
    "sampleSet",
    id
)
samples <- reactive({
    x <- sampleSet$assignments()
    req(x, is.data.table(x), nrow(x) > 0)   
    x 
})
assembleDataTable <- function(fileType, loadFn, colNames = NULL){
    x <- samples()
    do.call(rbind, lapply(seq_len(nrow(x)), function(i){
        dt <- loadFn(getSourceFilePath(x$Source_ID[i], fileType))
        if(!is.null(colNames)) setnames(dt, colNames)
        cbind(
            sourceId = x$Source_ID[i],
            sample = getSampleNames(sampleUniqueIds = paste(x$Project[i], x$Sample_ID[i], sep = ":")),
            dt
        )
    }))
}
getAmpliconKeys <- function(dt){
    paste(dt$sample, dt$amplicon, sep = ":")
}
amplicons <- reactive({
    x <- assembleDataTable("amplicons", fread, colNames = c(
        "amplicon","type","nReadPairs",
        "chrom1","side1","pos1","ref1","primer1",
        "chrom2","side2","pos2","ref2","primer2"
    ))
    x[, size := pos2 - pos1 + 1]
})
selectedAmplicons <- reactive({
    I <- ampliconsTable$rows_selected()
    req(I)
    amplicons()[I]    
})
moleculeTypes <- reactive({
    assembleDataTable("moleculeTypes", readRDS)
# [1] "sourceId"   "sample"     "molTypeId"  "amplicon"   "mergeLevel"  
# [6] "overlap"    "isRef"      "nReadPairs" "nMols"      "tLen"
# [11] "nMBases"    "nDBases"    "nIBases"    "avgQual"    "molClass"
# [16] "nAlns"      "nJxns"      "jxnsKey"    "maxIntMapQ" "alns"
# [21] "jxns"       "indels"     "seq1"       "qual1"      "seq2"
# [26] "qual2"      "jxnKeys"  
})
filteredMoleculeTypes <- reactive({
    moleculeTypes <- moleculeTypes()
    moleculeTypes[
        getAmpliconKeys(moleculeTypes) %in% getAmpliconKeys(selectedAmplicons()) & 
        mergeLevel > 0 & 
        grepl(CONSTANTS$svTypes$Deletion, jxnsKey, ignore.case = TRUE) # TODO: deploy dynamic SV type filtering
    ]
})
tableFilteredMoleculeTypes <- reactive({
    junctionI <- junctionsTable$rows_selected()
    if(!isTruthy(junctionI)) return(filteredMoleculeTypes())
    jxnKey <- filteredJunctions()[junctionI, jxnKey]
    filteredMoleculeTypes()[sapply(jxnKeys, function(x) jxnKey %in% x)]
})
junctions <- reactive({
    assembleDataTable("junctions", readRDS)
#  [1] "sourceId"   "sample"     "jxnKey"     "amplicon"   "molTypeIds"
#  [6] "nMolTypes"  "nMols"      "nReadPairs" "hasMerged"  "tLens"
# [11] "overlap"    "jxnBases"   "callTypes"  "svType"     "svSize"
# [16] "chrom1"     "side1"      "pos1"       "chrom2"     "side2"
# [21] "pos2"       "jxnId"
})
filteredJunctions <- reactive({
    junctions <- junctions()
    junctions[
        getAmpliconKeys(junctions) %in% getAmpliconKeys(selectedAmplicons()) &
        hasMerged == TRUE &
        svType == CONSTANTS$svTypes$Deletion &  # TODO: deploy dynamic SV type filtering
        TRUE
    ]
})

#----------------------------------------------------------------------
# construct the output tables
#----------------------------------------------------------------------
ampliconsTable <- bufferedTableServer(
    "ampliconsTable",
    id,
    input,
    tableData = reactive({
        amplicons()[, .SD, .SDcols = c(
            "sample",
            "amplicon","type",#"nReadPairs",
            "chrom1","side1","pos1",
            "chrom2","side2","pos2",
            "size"
        )]
    }),
    selection = 'multiple'
)
moleculeTypesTable <- bufferedTableServer(
    "moleculeTypesTable",
    id,
    input,
    tableData = reactive({
        tableFilteredMoleculeTypes()[, .SD, .SDcols = c(
            "sample",
            "amplicon","molTypeId",
            "nMols","nReadPairs",        
            "isRef","mergeLevel","overlap","tLen",
            "molClass","jxnsKey","avgQual","maxIntMapQ",
            "nAlns","nJxns"
        )]
    }),
    selection = 'single'
)
junctionsTable <- bufferedTableServer(
    "junctionsTable",
    id,
    input,
    tableData = reactive({
        filteredJunctions()[, .SD, .SDcols = c(
            "sample",
            "amplicon","nMolTypes","nMols","nReadPairs",
            "hasMerged","tLens","overlap","jxnBases","callTypes",
            "svType","svSize",
            "chrom1","side1","pos1",
            "chrom2","side2","pos2"
        )]
    }),
    selection = 'single'
)

#----------------------------------------------------------------------
# make junction position density plots
#----------------------------------------------------------------------
getDensity <- function(allValues, values){
    density <- merge(
        data.table(x = allValues),
        data.table(x = values)[, .(y = .N), by = "x"],
        by = "x",
        all.x = TRUE
    )
    density[is.na(y), y := 0]
    density[, .(
        x = x,
        y = y / sum (y)
    )]
}
plotData <- reactive({
    amplicons <- selectedAmplicons()
    filteredJunctions <- filteredJunctions()
    ampKeys <- getAmpliconKeys(amplicons)
    posLim <- c(min(amplicons$pos1), max(amplicons$pos2))
    sizeMax <- max(amplicons$pos2 - amplicons$pos1 + 1)
    list(
        amplicons = amplicons,
        ampKeys = ampKeys,
        triangle = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            filteredJunctions[getAmpliconKeys(filteredJunctions) == ampKeys[i], {
                size <- pos2 - pos1 + 1
                .(
                    x = pos1 + size / 2, # center of SV span
                    y = size             # SV size
                )
            }]
        }),
        positionDensity = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            filteredJunctions[
                getAmpliconKeys(filteredJunctions) == ampKeys[i], 
                getDensity(posLim[1]:posLim[2], c(pos1, pos2))
            ]
        }),
        sizeDensity = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            filteredJunctions[
                getAmpliconKeys(filteredJunctions) == ampKeys[i], 
                getDensity(0:sizeMax, pos2 - pos1 + 1)
            ]
        }),
        posLim = posLim, # in bp,
        sizeMax = sizeMax
    )
})
svTrianglePlot <- staticPlotBoxServer(
    "svTrianglePlot",
    maxHeight = "400px",
    points   = TRUE,
    legend  = TRUE,
    margins = FALSE,
    create = function() {
        d <- plotData()

        par(mar = c(4.1, 4.1, 0.1, 0.1))
        svTrianglePlot$initializeFrame(
            xlim = d$posLim / 1e6,
            ylim = c(0, d$sizeMax),
            xlab = "SV Center (Mb)",
            ylab = "SV Size (bp)"
        )

        for(i in seq_along(d$triangle)){
            svTrianglePlot$addPoints(
                x = d$triangle[[i]]$x / 1e6,
                y = d$triangle[[i]]$y,
                col = CONSTANTS$plotlyColors[[i]]
            )            
        }
        svTrianglePlot$addLegend(
            legend = d$ampKeys,
            col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
        )
    }
)
positionDensityPlot <- staticPlotBoxServer(
    "positionDensityPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = FALSE,
    create = function() {
        d <- plotData()

        par(mar = c(4.1, 4.1, 0.1, 0.1))
        positionDensityPlot$initializeFrame(
            xlim = d$posLim / 1e6,
            ylim = c(0, 0.05),
            xlab = "Endpoint Coordinate (Mb)",
            ylab = "Frequency"
        )

        for(i in seq_along(d$positionDensity)){
            svTrianglePlot$addLines(
                x = d$positionDensity[[i]]$x / 1e6,
                y = d$positionDensity[[i]]$y,
                col = CONSTANTS$plotlyColors[[i]]
            )            
        }
        svTrianglePlot$addLegend(
            legend = d$ampKeys,
            col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
        )
    }
)
sizeDensityPlot <- staticPlotBoxServer(
    "sizeDensityPlot",
    maxHeight = "400px",
    lines   = TRUE,
    legend  = TRUE,
    margins = FALSE,
    create = function() {
        d <- plotData()

        par(mar = c(4.1, 4.1, 0.1, 0.1))
        positionDensityPlot$initializeFrame(
            xlim = c(0, d$sizeMax),
            ylim = c(0, 0.05),
            xlab = "SV Size (bp)",
            ylab = "Density"
        )

        for(i in seq_along(d$sizeDensity)){
            svTrianglePlot$addLines(
                x = d$sizeDensity[[i]]$x,
                y = d$sizeDensity[[i]]$y,
                col = CONSTANTS$plotlyColors[[i]]
            )            
        }
        svTrianglePlot$addLegend(
            legend = d$ampKeys,
            col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
        )
    }
)

#----------------------------------------------------------------------
# make junction position density plots
#----------------------------------------------------------------------
getPremergedReads <- function(molType){
    x <- read_yaml(expandSourceFilePath(molType[, sourceId], "package.yml"))$task$collate$output
    file <- paste0(x[['data-name']], '.indexed_reads.bgz')
    file <- file.path(x[['output-dir']], x[['data-name']], file)
    if(!file.exists(file)) return(NULL) # this only works on an HPC server
    tbx <- Rsamtools::TabixFile(file)
    index <- floor(molType[, molTypeId] / 100) + 1 # molTypeId is one example of a read pair with the sequence/type
    i <- molType[, molTypeId] - (index - 1) * 100 + 1
    x <- Rsamtools::scanTabix(tbx, param = GenomicRanges::GRanges("X", IRanges::IRanges(index, width = 1)))
    x <- strsplit(x[[1]][i], "\t")[[1]]
    data.table(
        key_ = c("SEQ1","QUAL1","SEQ2","QUAL2"),
        value_ = x[4:7]
    )
}
moleculeTypeData <- reactive({
    I <- moleculeTypesTable$rows_selected()
    req(I)
    moleculeType <- tableFilteredMoleculeTypes()[I]
    dt <- data.table(
        key_ = names(moleculeType), 
        value_ = as.character(unlist(moleculeType))
    )    
    if(moleculeType[, mergeLevel > 0]){
        premerged <- getPremergedReads(moleculeType)
        if(!is.null(premerged)) dt <- rbind(dt, premerged)
    }
    dt   
})
moleculeTypeExpansion <- bufferedTableServer(
    "moleculeTypeExpansion",
    id,
    input,
    tableData = reactive({
        moleculeTypeData()
        # steps to creating an alignment of two merged reads
        # identify which read is R, which is L
        # rc(L read)
        # pad with spaces at right end of and left end of L
        # should always match the length of merged?
        # create a match line, | for match
        # create the output lines

        # steps to creating an alignment of reads to the amplicon
        # retrieve the amplicon sequence
        # parse the CIGAR strings for each alignment left to right
        # use to count number of times each query base is covered by an alignment
        # 0 = novel insertion, 1 = mapped, 2 = microhomology
        # or use Smith-Waterman (in genomex)?
        # quite possibly see svCapture, has some code (but not necessarily handled the same)
    }),
    selection = 'single',
    options = list(
        paging = FALSE
    )
)
moleculeDotPlot <- staticPlotBoxServer(
    "moleculeDotPlot",
    maxHeight = "400px",
    points   = TRUE,
    legend  = FALSE,
    margins = FALSE,
    create = function() {
        d <- moleculeTypeData()
        amplicon <- selectedAmplicons()[
            sourceId == d[key_ == "sourceId", value_] & 
            amplicon == d[key_ == "amplicon", value_]
        ]
        if(amplicon[, type != "expectOverlap"]) return(NULL)
        ref <- amplicon[, ref1]
        alns   <- lapply(d[key_ == "alns",   strsplit(value_, ":::")[[1]]], function(x) strsplit(x, ":")[[1]])
        jxns   <- lapply(d[key_ == "jxns",   strsplit(value_, ":::")[[1]]], function(x) strsplit(x, ":")[[1]])  # one jxn between every two alns
        indels <- lapply(d[key_ == "indels", strsplit(value_, ":::")[[1]]], function(x) strsplit(x, "::")[[1]]) # one entry per aln, with potential multiple jxns
        str(alns)

        par(mar = c(4.1, 4.1, 0.1, 0.1))
        positionDensityPlot$initializeFrame(
            xlim = c(amplicon$pos1, amplicon$pos2),
            ylim = c(1, as.integer(d[key_ == "tLen", value_])),
            xlab = "Reference Coordinate",
            ylab = "Position in Query Molecule"
        )        


        return(NULL)

        # alns = 0:chr1:33764688:60:331M-48D-196M
        #        2048:chr1:33764688:60:128M447S:::0:chr1:33764951:60:263S312M
        # jxns   1:D:32:chr1/L/33765192:chr1/R/33765225:-32:GATAATTTTTTTTTTTCTTTTTTTTTTAGTTGGAAAT
        # indels 1:D:48:chr1/L/33765018:chr1/R/33765067:NA:NA





        for(i in seq_along(d$sizeDensity)){
            svTrianglePlot$addLines(
                x = d$sizeDensity[[i]]$x,
                y = d$sizeDensity[[i]]$y,
                col = CONSTANTS$plotlyColors[[i]]
            )            
        }
        svTrianglePlot$addLegend(
            legend = d$ampKeys,
            col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
        )
    }
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
        # outcomes <<- listToReactiveValues(bm$outcomes)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = reactive({ list(
    ) }),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
