#----------------------------------------------------------------------
# handle junction cluster loading and filtering
#----------------------------------------------------------------------

# parse a track's item's into sourceId -> samples list
getSourceIdFromSample_svPore <- function(sample){
    uploadName <- appStepNamesByType$upload
    samples <- as.data.table(app[[uploadName]]$outcomes$samples())        
    samples[Sample_ID == sample$Sample_ID & Project == sample$Project, Source_ID]
}
getSvPoreSampleSources <- function(samples){ # samples is a list of sample lists, with Sample_ID and Project
    sources <- list()
    for(sample in samples){
        sourceId <- getSourceIdFromSample_svPore(sample)
        sources[[sourceId]] <- rbind(sources[[sourceId]], as.data.table(sample))
    }    
    sources
}

# load all junction clusters for a specific sourceId (not filtered by region, type, or sample yet)
loadJunctionClusters <- function(sourceId){
    req(sourceId)
    svPoreCache$get(
        'junctionClusters', 
        key = sourceId, 
        permanent = TRUE, 
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "loading source JCs")
            jc <- readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
            cm <- readRDS(getSourceFilePath(sourceId, "chromosomesFile"))
            centerPos <- function(p1, p2) pmin(p1, p2) + abs(p2 - p1) / 2
            jc[, ":="(
                cChrom1 = unlist(cm$revChromIndex[cChromIndex1]),
                cChrom2 = unlist(cm$revChromIndex[cChromIndex2]),
                nodeCenter   = centerPos(abs(node1), abs(node2)),
                refPosCenter = centerPos(cRefPos1,   cRefPos2),
                size = abs(abs(node2) - abs(node1)), # eventSize has the chrom-level size, where translocation size == 0
                samples = paste0(",", samples, ",") # for simpler matching of single samples
            )]
            jc
        }
    )$value
}

# filter junction clusters by sample, e.g., specified by browser items
filterJCsBySample <- function(sourceId, samples_){
    svPoreCache$get(
        'junctionClusters', 
        keyObject = list(sourceId = sourceId, samples = samples_), 
        permanent = TRUE, 
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            jc <- loadJunctionClusters(sourceId)
            startSpinner(session, message = "filtering JCs by sample")
            I <- FALSE
            samples_ <- paste0(",", samples_, ",")
            for(sample_ in samples_) I <- I | jc[, grepl(sample_, samples)]
            jc[I]
        }
    )$value
}

# filter junction clusters based on track or other settings
svPoreFilterDefaults <- list(
    Min_SV_Size = 1,
    Max_SV_Size = 0,
    Min_Samples_With_SV = 1,
    Max_Samples_With_SV = 0,
    Min_Source_Molecules = 1,
    Max_Source_Molecules = 0,    
    SV_Type = c("Del","Dup","Inv"),
    Show_ChrM = "never"
)
applySettingsToJCs <- function(sourceId, samples, track){
    svPoreCache$get(
        'junctionClusters', 
        keyObject = list(sourceId = sourceId, samples = samples, settings = track$settings$all()), 
        permanent = FALSE,
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            jc <- filterJCsBySample(sourceId, samples)
            startSpinner(session, message = paste("applying JC settings"))
            filters <- track$settings$Filters()
            filters <- lapply(names(svPoreFilterDefaults), function(filter){
                if(is.null(filters[[filter]])) svPoreFilterDefaults[[filter]]
                else if(!is.null(filters[[filter]]$selected)) filters[[filter]]$selected else filters[[filter]]$value
            })
            names(filters) <- names(svPoreFilterDefaults)

            if(filters$Min_SV_Size > 1) jc <- jc[size >= filters$Min_SV_Size]
            if(filters$Max_SV_Size > 0) jc <- jc[size <= filters$Max_SV_Size] #  & edgeType != "T"

            if(filters$Min_Samples_With_SV > 1) jc <- jc[nSamples >= filters$Min_Samples_With_SV]
            if(filters$Max_Samples_With_SV > 0) jc <- jc[nSamples <= filters$Max_Samples_With_SV]

            if(filters$Min_Source_Molecules > 1) jc <- jc[nInstances >= filters$Min_Source_Molecules]
            if(filters$Max_Source_Molecules > 0) jc <- jc[nInstances <= filters$Max_Source_Molecules]

            if(length(filters$SV_Type) > 0) jc <- jc[svPore$jxnTypes[edgeType, name] %in% filters$SV_Type]

            if(filters$Show_ChrM != "always"){
                hasChrM <- jc[, cChrom1 == "chrM" | cChrom2 == "chrM"]
                if(filters$Show_ChrM == "translocations_only") jc <- jc[hasChrM == FALSE | xor(cChrom1 == "chrM", cChrom2 == "chrM")]
                else if(filters$Show_ChrM == "never")          jc <- jc[hasChrM == FALSE]
            }

            jc %>% svPore$setJunctionPointColors(track) %>% svPore$setJunctionPointSizes(track)
        }
    )$value
}

# filter junction clusters by browser coordinates
filterJCsByRange <- function(jc, coord, rangeType, chromOnly = TRUE){
    startSpinner(session, message = "filtering junction clusters to window")
    isWholeGenome <- coord$chromosome == "all"
    if(!isWholeGenome && chromOnly) jc <- jc[
        cChrom1 == coord$chromosome & 
        cChrom2 == coord$chromosome
    ]
    switch(
        rangeType,
        center = {
            jc <- jc[between(if(isWholeGenome) as.numeric(nodeCenter) else refPosCenter, coord$start, coord$end)]
            jc[ , center := if(isWholeGenome) nodeCenter else refPosCenter]
            jc
        },
        endpoint = {
            jc[ , ":="(
                pos1 = if(isWholeGenome) abs(node1) else cRefPos1,
                pos2 = if(isWholeGenome) abs(node2) else cRefPos2
            )]
            jc[ , ":="(
                pos1In = (isWholeGenome | cChrom1 == coord$chrom) & between(as.numeric(pos1), coord$start, coord$end),
                pos2In = (isWholeGenome | cChrom2 == coord$chrom) & between(as.numeric(pos2), coord$start, coord$end)
            )]
            jc[pos1In | pos2In]
        },
        jc
    )   
}

# get a single junction cluster for the object table
getJunctionCluster <- function(x){
    loadJunctionClusters(x$sourceId)[clusterN == x$clusterN]
}
