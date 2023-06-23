#----------------------------------------------------------------------
# handle junction cluster loading and filtering
#----------------------------------------------------------------------

# load all junction clusters for a specific sourceId (not filtered by region, type, or sample yet)
loadJunctionClusters <- function(sourceId = NULL, sample = NULL){
    if(is.null(sourceId)) {
        uploadName <- appStepNamesByType$upload
        samples <- as.data.table(app[[uploadName]]$outcomes$samples())        
        sourceId <- samples[Sample_ID == sample$Sample_ID & Project == sample$Project, Source_ID]
    }
    req(sourceId)
    jc <- svPoreCache$get('junctionClusters', key = sourceId, permanent = TRUE, from = "disk", create = "always", createFn = function(...) {
        # startSpinner(session, message = "loading junction clusters")
        jc <- readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
        cm <- readRDS(getSourceFilePath(sourceId, "chromosomesFile"))
        centerPos <- function(p1, p2) pmin(p1, p2) + abs(p2 - p1) / 2
        jc[, ":="(
            cChrom1 = unlist(cm$revChromIndex[cChromIndex1]),
            cChrom2 = unlist(cm$revChromIndex[cChromIndex2]),
            nodeCenter   = centerPos(abs(node1), abs(node2)),
            refPosCenter = centerPos(cRefPos1,   cRefPos2),
            size = abs(abs(node2) - abs(node1)),
            nSampleInstances = nInstances # unless overridden by a sample filter, below
        )]
        # stopSpinner(session)
        jc
    })$value
}

# filter junction clusters by sample, e.g., specified by browser items
filterJCsBySample <- function(jc, sample){
    nSampleInstances_ <- jc[[sample$Sample_ID]]
    jc[, nSampleInstances := nSampleInstances_]
    jc[nSampleInstances > 0]
}

# filter junction clusters by browser coordinates
filterJCsByRange <- function(jc, coord, rangeType){
    isWholeGenome <- coord$chromosome == "all"
    if(!isWholeGenome) jc <- jc[
        cChrom1 == coord$chromosome & 
        cChrom2 == coord$chromosome
    ]
    jc <- switch(
        rangeType,
        center = jc[between(if(isWholeGenome) as.numeric(nodeCenter) else refPosCenter, coord$start, coord$end)],
        jc
    )
    jc    
}

# filter junction clusters based on track or other settings
filterJCsBySettings <- function(jc, settings){
    filters <- settings$SV_Filters()

    if(filters$Min_SV_Size$value > 1) jc <- jc[size >= filters$Min_SV_Size$value]
    if(filters$Max_SV_Size$value > 0) jc <- jc[size <= filters$Max_SV_Size$value]

    if(filters$Min_Samples_With_SV$value > 1) jc <- jc[nSamples >= filters$Min_Samples_With_SV$value]
    if(filters$Max_Samples_With_SV$value > 0) jc <- jc[nSamples <= filters$Max_Samples_With_SV$value]

    if(filters$Min_Source_Molecules$value > 1) jc <- jc[nInstances >= filters$Min_Source_Molecules$value]
    if(filters$Max_Source_Molecules$value > 0) jc <- jc[nInstances <= filters$Max_Source_Molecules$value]

    if(length(filters$SV_Type$value) > 1) jc <- jc[svPore$jxnTypes[edgeType, name] %in% filters$SV_Type$value]

    jc    
}
