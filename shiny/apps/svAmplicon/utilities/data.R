# support for svAmplicon data loading

# samples from a selected sampleSet
samplesReactive <- function(sampleSet) reactive({
    x <- sampleSet$assignments()
    req(x, is.data.table(x), nrow(x) > 0)   
    x 
})

# general function for loading a data.table from the pipeline output
# prepend columns sourceId and sample to the pipeline data.table
assembleDataTable <- function(samples, fileType, loadFn, colNames = NULL){
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

# amplicons; could be more than one for a given sample if a mixed library
getAmpliconKeys <- function(dt){
    req(dt, nrow(dt) > 0)
    ampId <- if("amplicon" %in% names(dt)) "amplicon" else "ampliconId"
    paste(dt$sample, dt[[ampId]], sep = ":")
}
ampliconsReactive <- function(samples) reactive({
    x <- assembleDataTable(samples, "amplicons", fread, colNames = c(
        "amplicon","type","nReadPairs",
        "chrom1","side1","pos1","ref1","primer1",
        "chrom2","side2","pos2","ref2","primer2"
    ))
    x[, size := pos2 - pos1 + 1]
})
ampliconsTableServer <- function(parentId, input, amplicons, selection = "single") bufferedTableServer(
    "ampliconsTable",
    parentId,
    input,
    tableData = reactive({
        amplicons()[, .SD, .SDcols = c(
            "sample",
            # "amplicon","type",#"nReadPairs",
            "chrom1","side1","pos1",
            "chrom2","side2","pos2",
            "size"
        )]
    }),
    selection = selection
)
selectedAmpliconsReactive <- function(amplicons, ampliconsTable) reactive({
    I <- ampliconsTable$rows_selected()
    req(I)
    amplicons()[I]    
})

# moleculeTypes, i.e., one row for each distinct path through an amplicon molecule
moleculeTypesReactive <- function(samples, selectedAmplicons) reactive({ 
    amplicons <- selectedAmplicons()
    req(amplicons)
    x <- assembleDataTable(samples, "moleculeTypes", readRDS) 
    x[isRef == TRUE, pathClass := "-"] # avoids having a blank list name
    x[, ":="(
        ampliconKey = getAmpliconKeys(x),
        sampleMolTypeKey = paste(sample, molKey, sep = "::")
    )]
    x[ampliconKey %in% getAmpliconKeys(amplicons)]
})

# pathClasses, i.e., one row for each kind of path (D, TT, VV, DID, etc.) through an amplicon molecule
pathClassesReactive <- function(moleculeTypes) reactive({
    x  <- moleculeTypes()
    pc <- CONSTANTS$groupedPathClasses
    x[, pathClass_ := if(pathClass %in% names(pc)) pc[[pathClass]] else "other", by = "pathClass"]    
    x[, .(
        nReadPairs = sum(nReadPairs),
        nMols = sum(nMols),         
        nMolTypes = .N
    ), by = .(pathClass_)]
})
pathClassesTableServer <- function(parentId, input, pathClasses, selection = "single") bufferedTableServer(
    "pathClassesTable",
    parentId,
    input,
    tableData = reactive({ pathClasses() }),
    options = list(
        paging = FALSE,
        searching = FALSE
    ),
    selection = selection
)
pathClassMoleculeTypesReactive <- function(moleculeTypes, pathClasses, pathClassesTable) reactive({
    I <- pathClassesTable$rows_selected()
    x <- if(isTruthy(I)) moleculeTypes()[pathClass_ == pathClasses()[I, pathClass_]]
         else moleculeTypes()  
    thresholds <- app$keepReject$thresholds()[[1]]         
    x[, passedQualityFilters := minBaseQual >= thresholds$minBaseQual & minMapQ >= thresholds$minMapQ]
    x
})
keepRejectMoleculeTypesReactive <- function(pathClassMoleculeTypes, input) reactive({
    x <- pathClassMoleculeTypes()
    req(x, input$moleculeTypeFilter)
    filterState <- if(input$moleculeTypeFilter == "Kept") TRUE else FALSE    
    x[passedQualityFilters == filterState]
})

# moleculeType, i.e., metadata on a single selected moleculeType
moleculeMetadataReactive <- function(selectedAmplicon, moleculeTypeStepper) reactive({
    amplicon <- selectedAmplicon()
    moleculeType <- moleculeTypeStepper$getCurrent()
    req(amplicon, moleculeType)

    # molecule-level information
    alnReadNs <- as.integer(strsplit(moleculeType$readNs, ":")[[1]])
    jxnReadNs <- integer()
    if(!moleculeType$isRef){
        for(i in 1:(length(alnReadNs) - 1)) if(alnReadNs[i] == alnReadNs[i + 1]) jxnReadNs <- c(jxnReadNs, alnReadNs[i])
    }
    moleculeType$tLen1 <- nchar(moleculeType$seq1)
    moleculeType$tLen2 <- if(moleculeType$mergeLevel > 0) NA else nchar(moleculeType$seq2)
    moleculeType$alnReadNs <- paste(alnReadNs, collapse = ":")
    moleculeType$jxnReadNs <- paste(jxnReadNs, collapse = ":")

    # out-of-amplicon segments
    chroms <- strsplit(moleculeType$chroms, ":")[[1]]
    poss   <- as.integer(strsplit(moleculeType$poss, ":")[[1]])
    ampliconic <- chroms == amplicon$chrom1 & poss >= amplicon$pos1 & poss <= amplicon$pos2

    list(
        moleculeType = moleculeType,
        amplicon = amplicon,
        ampliconSize = amplicon$pos2 - amplicon$pos1 + 1,
        tLen = moleculeType[,
            if(mergeLevel == 0) nchar(seq1) + nchar(seq2) - (if(is.na(overlap)) 0 else overlap) # sometimes overlap is NA
            else nchar(seq1)
        ],
        isMerged = moleculeType$mergeLevel > 0,
        ampliconic = ampliconic,

        # alignment-level information
        chroms = chroms,
        strands = strsplit(moleculeType$strands, ":")[[1]],
        poss   = poss,
        cigars = strsplit(moleculeType$cigars, ":")[[1]],
        mapQs = strsplit(moleculeType$mapQs, ":")[[1]],
        baseQuals = strsplit(moleculeType$baseQuals, ":")[[1]],
        alnReadNs = alnReadNs,
        jxnReadNs = jxnReadNs,
        sizes  = as.integer(strsplit(moleculeType$sizes, ":")[[1]]),
        insSizes = as.integer(strsplit(moleculeType$insSizes, ":")[[1]]),
        types = strsplit(moleculeType$pathClass, "")[[1]]
    )
})
moleculeMetadataUI <- function(moleculeMetadata) renderUI({
    mmd <- moleculeMetadata()
    req(mmd)
    mt <- mmd$moleculeType
    keyValuePair <- function(key) {
        tags$div(
            tags$div(paste(key, ":"), style = "display: inline-block; width: 100px; text-align: right; margin-right: 5px;"),
            tags$div(gsub(":", " ", mt[[key]]), style = "display: inline-block; margin-left:  5px;")
        )
    }
    tags$div(
        lapply(c(
            "ampliconId","molId",
            "nReadPairs","nMols",
            "mergeLevel","overlap","tLen1","tLen2",
            "pathClass",
            "alnReadNs",
            "strands","chroms","poss","cigars","mapQs","baseQuals",
            "jxnReadNs",
            "sizes","insSizes"
        ), keyValuePair)
    )
})

# junctions, i.e., one row for each distinct SV edge in an amplicon molecule
junctionsReactive <- function(samples, moleculeTypes, selectedAmplicons) reactive({ 
    amplicons <- selectedAmplicons()
    moleculeTypes <- moleculeTypes()
    req(amplicons, moleculeTypes)
    x <- assembleDataTable(samples, "junctions", readRDS) 
    x <- x[edgeClass != CONSTANTS$edgeTypes$ALIGNMENT]
    x[, jxnKey := paste(sample, ampliconId, nodePair, insSize, sep = "::")]
    expandedJunctions <- x[, .(sampleMolTypeKey = paste(sample, unlist(molTypeKeys), sep = "::")), by = .(jxnKey)]
    thresholds <- app$keepReject$thresholds(amplicons)
    moleculeTypes <- do.call(rbind, lapply(getAmpliconKeys(amplicons), function(aKey){
        ts <- thresholds[[aKey]]
        moleculeTypes[ampliconKey == aKey & minBaseQual >= ts$minBaseQual & minMapQ >= ts$minMapQ]
    }))
    jxnKeys <- expandedJunctions[sampleMolTypeKey %in% moleculeTypes$sampleMolTypeKey, unique(jxnKey)]
    x <- x[jxnKey %in% jxnKeys] 
    jt <- CONSTANTS$junctionTypes    
    x[, junctionType := if(edgeClass %in% names(jt)) jt[[edgeClass]] else "other", by = "edgeClass"]    
    x
})

# junction Types, i.e., one row for each kind of junction (D, T, V, etc.)
junctionTypesReactive <- function(junctions) reactive({
    junctions <- junctions()
    req(junctions)
    junctions[, .(
        nReadPairs = sum(nReadPairs),
        nMols = sum(nMol),         
        nMolTypes = sum(nMolTypes),
        nJunctions = .N
    ), by = .(junctionType)]
})
junctionTypesTableServer <- function(parentId, input, junctionTypes, selection = "multiple") bufferedTableServer(
    "junctionsTypesTable",
    parentId,
    input,
    tableData = reactive({ junctionTypes() }),
    options = list(
        paging = FALSE,
        searching = FALSE
    ),
    selection = selection
)
junctionTypesJunctionsReactive <- function(junctions, junctionTypes, junctionTypesTable) reactive({
    I <- junctionTypesTable$rows_selected()
    x <- if(isTruthy(I)) junctions()[junctionType %in% junctionTypes()[I, junctionType]]
         else junctions()  
    x
})
junctionsTableServer <- function(parentId, input, junctionPlotData) bufferedTableServer(
    "junctionsTable",
    parentId,
    input,
    tableData = reactive({ junctionPlotData()$junctions[, .SD, .SDcols = c(
        "sample",
        "ampliconId",
        # "nodePair",
        "chrom1",
        "strand1",
        "pos1",
        "chrom2",
        "strand2",
        "pos2",
        "insSize", 
        "nReadPairs", 
        "nMol", 
        "nMolTypes", 
        "junctionType", 
        "size"
    )] }),
    selection = "none"
)

    # # ampliconId
    # # nodePair
    # # insSize
    # nReadPairs = sum(nReadPairs),
    # nMol = length(unique(molKey)), # beware of possibility of rare molecules with a recurring junction
    # nMolTypes = length(unique(molKey[isIndexMol])), 
    # edgeClass = edgeClass[1],
    # mergeLevel = max(mergeLevel),
    # size = size[1],
    # mapQ = max(mapQ),
    # molTypeKeys = list(unique(molKey[isIndexMol]))


