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
        sampleUniqueId <- paste(x$Project[i], x$Sample_ID[i], sep = ":") # thus, sampleUniqueId = Project:Sample_ID
        cbind(
            sourceId       = x$Source_ID[i],
            sampleUniqueId = sampleUniqueId, # for tracking
            sampleName     = getSampleNames(sampleUniqueIds = sampleUniqueId), # for display
            dt
        )
    }))
}

# amplicons; could be more than one for a given sample if a mixed library
getAmpliconKeys <- function(dt){
    req(dt, nrow(dt) > 0)
    ampId <- if("amplicon" %in% names(dt)) "amplicon" else "ampliconId"
    paste(dt$sampleUniqueId, dt[[ampId]], sep = ":") # thus, ampliconKey = Project:Sample_ID:amplicon = sampleUniqueId:amplicon
}
CONSTANTS$ampliconsColNames <- c(
    "amplicon","type","nReadPairs",
    "chrom1","side1","pos1","ref1","primer1",
    "chrom2","side2","pos2","ref2","primer2"
)
fillAmpliconDerivedColumns <- function(amplicons){
    amplicons[, size := if(chrom1 == chrom2) pos2 - pos1 + 1 else NA_integer_]
    if("side1" %in% names(amplicons)){
        amplicons[, ":="(
            strand1 = if(side1 == "R") "+" else "-",
            strand2 = if(side2 == "R") "+" else "-"
        )]
    } else {
        amplicons[, ":="(
            side1 = if(strand1 == "+") "R" else "L",
            side2 = if(strand2 == "+") "R" else "L"    
        )]
    }
    amplicons
}
ampliconsReactive <- function(samples) reactive({
    assembleDataTable(samples, "amplicons", fread, colNames = CONSTANTS$ampliconsColNames) %>%
    fillAmpliconDerivedColumns()
})
ampliconsTableServer <- function(parentId, input, amplicons, selection = "single") bufferedTableServer(
    "ampliconsTable",
    parentId,
    input,
    tableData = reactive({
        amplicons()[, .SD, .SDcols = c(
            "sampleName",
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
    startSpinner(session, message = "loading molecule types")
    x <- assembleDataTable(samples, "moleculeTypes", readRDS) 
    x[isReference == TRUE, pathEdgeTypes := "-"] # avoids having a blank list name
    x[, ":="(
        ampliconKey = getAmpliconKeys(x),
        sampleMolTypeKey = paste(sampleUniqueId, molKey, sep = ":") # thus, sampleMolTypeKey = Project:Sample_ID:amplicon:indexMolId = ampliconKey:indexMolId
    )]
    x <- x[ampliconKey %in% getAmpliconKeys(amplicons)]
    stopSpinner(session)
    x
})

# pathClasses, i.e., one row for each kind of edgeTypes path (D, TT, VV, DID, etc.) through an amplicon molecule
pathClassesReactive <- function(moleculeTypes) reactive({
    x  <- moleculeTypes()
    startSpinner(session, message = "loading path classes")
    pc <- CONSTANTS$groupedPathClasses
    x[, pathClass := if(pathEdgeTypes %in% names(pc)) pc[[pathEdgeTypes]] else "other", by = "pathEdgeTypes"]    
    x  <- x[, .(
        nReadPairs = sum(nReadPairs),
        nMolecules = sum(nMolecules),         
        nMolTypes = .N
    ), by = .(pathClass)]
    stopSpinner(session)
    x
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
pathClassMoleculeTypesReactive <- function(moleculeTypes, pathClasses, pathClassesTable, selectedAmplicons) reactive({
    I <- pathClassesTable$rows_selected()
    x <- if(isTruthy(I)) moleculeTypes()[pathClass %in% pathClasses()[I, pathClass]]
         else moleculeTypes()    
    thresholds <- app$keepReject$thresholds(selectedAmplicons) 
    for(aKey in names(thresholds)){
        ts <- thresholds[[aKey]]
        x[ampliconKey == aKey, passedQualityFilters := minBaseQual >= ts$minBaseQual & minMapQ >= ts$minMapQ]        
    }
    x
})
keepRejectMoleculeTypesReactive <- function(moleculeTypes, input, selectedAmplicon, kept) reactive({
    x <- moleculeTypes()
    req(x, input$moleculeTypeFilter)
    filterState <- if(input$moleculeTypeFilter == "Kept") TRUE else FALSE  
    ampliconKey <- getAmpliconKeys(selectedAmplicon())
    userRejects <- if(is.null(kept[[ampliconKey]])) character() else na.omit(sapply(names(kept[[ampliconKey]]), function(sampleMolTypeKey){
        if(kept[[ampliconKey]][[sampleMolTypeKey]]) NA else sampleMolTypeKey
    }))
    if(filterState == TRUE) x[passedQualityFilters == TRUE & !(sampleMolTypeKey %in% userRejects)]
                       else x[passedQualityFilters == FALSE |  sampleMolTypeKey %in% userRejects]
})

# moleculeType, i.e., metadata on a single selected moleculeType
moleculeMetadataReactive <- function(selectedAmplicon, moleculeTypeStepper) reactive({
    amplicon <- selectedAmplicon()
    moleculeType <- moleculeTypeStepper$getCurrent()
    req(amplicon, moleculeType)

    # molecule-level information
    moleculeType$tLen1 <- nchar(moleculeType$seq1)
    moleculeType$tLen2 <- if(moleculeType$mergeLevel > 0) NA else nchar(moleculeType$seq2)

    # out-of-amplicon segments
    split_ <- function(column, isInteger = FALSE){
        x <- strsplit(moleculeType[[column]], ":")[[1]]
        if(isInteger) x <- as.integer(x)
        x
    }
    cigars <- split_("cigars")
    nAlns <- length(cigars)
    nJxns <- nAlns + 1
    chroms <- split_("chroms")
    poss   <- split_("poss", TRUE)
    alnPos <- poss[seq(1, length(poss), 2)]
    padding <- 10
    ampliconic <- if(nAlns <= 2) rep(TRUE, nAlns) 
                  else chroms == amplicon$chrom1 & alnPos >= amplicon$pos1 - padding & alnPos <= amplicon$pos2 + padding
    list(
        # molecule-level information
        moleculeType = moleculeType, # unparsed information, some repeated below for convenient access
        amplicon = amplicon,
        ampliconSize = amplicon$pos2 - amplicon$pos1 + 1,
        tLen = moleculeType[,
            if(mergeLevel == 0) nchar(seq1) + nchar(seq2) - (if(is.na(overlap)) 0 else overlap) # sometimes overlap is NA
            else nchar(seq1)
        ],
        isMerged = moleculeType$mergeLevel > 0,

        # node-level information
        poss        = poss,        

        # junction-level information
        edgeTypes   = split_("edgeTypes"), # fields that have the same value for the entire junction
        eventSizes  = split_("eventSizes", TRUE),
        insertSizes = split_("insertSizes", TRUE),
        jxnBases    = split_("jxnBases"),
        insertBases = split_("insertBases"),
        mapQs       = split_("mapQs", TRUE),
        baseQuals   = split_("baseQuals", TRUE),

        # alignment-level information, one value per alignment
        chroms      = chroms,
        strands     = split_("strands"),
        cigars      = cigars,
        readNs      = split_("readNs", TRUE),         
        ampliconic = ampliconic       
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
            "nReadPairs","nMolecules",
            "mergeLevel","overlap","tLen1","tLen2",
            "poss",
            "chroms","strands","readNs","cigars",
            "edgeTypes","eventSizes","insertSizes",
            "mapQs","baseQuals"
        ), keyValuePair)
    )
})

# junctions, i.e., one row for each distinct SV edge in an amplicon molecule
junctionsReactive <- function(samples, selectedAmplicons, pathClassMoleculeTypes) reactive({ 
    amplicons <- selectedAmplicons()
    moleculeTypes <- pathClassMoleculeTypes()
    req(amplicons, moleculeTypes)
    startSpinner(session, message = "loading junctions")
    junctions <- assembleDataTable(samples, "junctions", readRDS) 
    junctions[, jxnKey := paste(sampleUniqueId, ampliconId, node1, node2, insertSize, insertBases, sep = ":")]
    expandedJunctions <- junctions[, .(sampleMolTypeKey = paste(sampleUniqueId, unlist(molTypeKeys), sep = ":")), by = .(jxnKey)]
    thresholds <- app$keepReject$thresholds(amplicons)
    kept <- app$keepReject$outcomes()$kept
    moleculeTypes <- do.call(rbind, lapply(getAmpliconKeys(amplicons), function(aKey){
        ts <- thresholds[[aKey]]
        userRejects <- if(is.null(kept[[aKey]])) character() else na.omit(sapply(names(kept[[aKey]]), function(sampleMolTypeKey){
            if(kept[[aKey]][[sampleMolTypeKey]]) NA else sampleMolTypeKey
        }))
        moleculeTypes[
            ampliconKey == aKey & 
            minBaseQual >= ts$minBaseQual & 
            minMapQ >= ts$minMapQ & 
            !(sampleMolTypeKey %in% userRejects)
        ]
    }))
    jxnKeys <- expandedJunctions[sampleMolTypeKey %in% moleculeTypes$sampleMolTypeKey, unique(jxnKey)]
    junctions <- junctions[jxnKey %in% jxnKeys] 
    jt <- CONSTANTS$junctionTypes    
    junctions[, jxnType := if(edgeType %in% names(jt)) jt[[edgeType]] else "other", by = "edgeType"]  
    stopSpinner(session)  
    junctions
})

# junction Types, i.e., one row for each kind of junction (D, T, V, etc.)
junctionTypesReactive <- function(junctions) reactive({
    junctions <- junctions()
    req(junctions)
    junctions[, .(
        nReadPairs = sum(nReadPairs),
        nMolecules = sum(nMolecules),         
        nMolTypes = sum(nMolTypes),
        nJunctions = .N
    ), by = .(jxnType)]
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
    x <- if(isTruthy(I)) junctions()[jxnType %in% junctionTypes()[I, jxnType]]
         else junctions()  
    x
})
junctionsTableServer <- function(parentId, input, junctionPlotData) bufferedTableServer(
    "junctionsTable",
    parentId,
    input,
    tableData = reactive({ junctionPlotData()$junctions[, .SD, .SDcols = c(
        "sampleName",
        "ampliconId",
        "jxnType",         
        "chrom1",
        "strand1",
        "pos1",
        "chrom2",
        "strand2",
        "pos2",
        "eventSize",
        "insertSize",
        "jxnBases",
        "nReadPairs", 
        "nMolecules", 
        "nMolTypes"
    )] }),
    selection = "none"
)
