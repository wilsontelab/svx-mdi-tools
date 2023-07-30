#----------------------------------------------------------------------
# handle unique junction loading and filtering
# these are generic functions that apply to potentially any SVX app
# expects:
#   sessionCache
#   chromosomesFile provided in data.package
#   loadFn(sourceId) that returns a data.table with appropriately named display and filter columns
#----------------------------------------------------------------------

# load all unique junctions for a specific sourceId (not filtered by region, type, or sample yet)
svx_loadJunctions <- function(sourceId, loadFn){
    req(sourceId)
    sessionCache$get(
        'junctions', 
        key = sourceId, 
        permanent = TRUE, 
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "loading junctions")
            jxns <- loadFn(sourceId)
            chroms <- readRDS(getSourceFilePath(sourceId, "chromosomesFile")) 
            getCenter <- function(p1, p2) pmin(p1, p2) + abs(p2 - p1) / 2  
            jxns[, ":="(
                cChrom1 = unlist(chroms$revChromIndex[cChromIndex1]),
                cChrom2 = unlist(chroms$revChromIndex[cChromIndex2]),
                nodeCenter   = getCenter(abs(node1), abs(node2)),
                refPosCenter = getCenter(cRefPos1,   cRefPos2),
                size = abs(abs(node2) - abs(node1)), # eventSize has the chrom-level size, where translocation size == 0
                samples = paste0(",", samples, ",") # for simpler matching of single samples
            )]    
            jxns
        }
    )$value
}

# filter unique junctions by sample, e.g., specified by browser items
svx_filterJunctionsBySample <- function(sourceId, samples_, loadFn){
    sessionCache$get(
        'junctions', 
        keyObject = list(sourceId = sourceId, samples = samples_), 
        permanent = TRUE, 
        from = "ram",
        create = "asNeeded", 
        createFn = function(...) {
            jxns <- svx_loadJunctions(sourceId, loadFn)
            startSpinner(session, message = "filtering junctions")
            I <- FALSE
            samples_ <- paste0(",", samples_, ",")
            for(sample_ in samples_) I <- I | jxns[, grepl(sample_, samples)]
            jxns[I]
        }
    )$value
}

# filter junction clusters based on track or other settings
svx_filterJunctionsBySettings <- function(track, sourceId, samples, loadFn){
    jxns <- sessionCache$get(
        'junctions', 
        keyObject = list(sourceId = sourceId, samples = samples, settings = track$settings$all()), 
        permanent = FALSE,
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            jxns <- svx_filterJunctionsBySample(sourceId, samples, loadFn)
            startSpinner(session, message = paste("getting junctions"))

            filters <- track$settings$Filters()
            filters <- lapply(names(svx_filterDefaults), function(filter){
                if(is.null(filters[[filter]])) svx_filterDefaults[[filter]]
                else if(!is.null(filters[[filter]]$selected)) filters[[filter]]$selected else filters[[filter]]$value
            })
            names(filters) <- names(svx_filterDefaults)

            if(filters$Min_SV_Size > 1) jxns <- jxns[size >= filters$Min_SV_Size]
            if(filters$Max_SV_Size > 0) jxns <- jxns[size <= filters$Max_SV_Size] #  & edgeType != svx_edgeTypes$TRANSLOCATION

            jxns <- jxns[insertSize >= filters$Min_Insert_Size] # can be negative
            jxns <- jxns[insertSize <= filters$Max_Insert_Size]

            if(filters$Min_Samples_With_SV > 1) jxns <- jxns[nSamples >= filters$Min_Samples_With_SV]
            if(filters$Max_Samples_With_SV > 0) jxns <- jxns[nSamples <= filters$Max_Samples_With_SV]

            if(filters$Min_Source_Molecules > 1) jxns <- jxns[nInstances >= filters$Min_Source_Molecules]
            if(filters$Max_Source_Molecules > 0) jxns <- jxns[nInstances <= filters$Max_Source_Molecules]

            if(filters$Min_Sequenced_Molecules > 0) jxns <- jxns[nSequenced >= filters$Min_Sequenced_Molecules]
            if(filters$Max_Linked_Junctions > 0) jxns <- jxns[nLinkedJunctions <= filters$Max_Linked_Junctions]

            if(filters$Min_Map_Quality > 0) jxns <- jxns[mapQ >= filters$Min_Map_Quality] 
            if(filters$Min_Flank_Length > 0) jxns <- jxns[flankLength >= filters$Min_Flank_Length] 

            if(length(filters$SV_Type) > 0) {
                edgeTypes <- svx_jxnType_nameToX(filters$SV_Type, "code")
                jxns <- jxns[edgeType %in% edgeTypes]
            }

            if(filters$Show_ChrM != "always"){
                hasChrM <- jxns[, cChrom1 == "chrM" | cChrom2 == "chrM"]
                if(filters$Show_ChrM == "translocations_only") jxns <- jxns[hasChrM == FALSE | xor(cChrom1 == "chrM", cChrom2 == "chrM")]
                else if(filters$Show_ChrM == "never")          jxns <- jxns[hasChrM == FALSE]
            }

            jxns %>% svx_setJunctionPointColors(track) %>% svx_setJunctionPointSizes(track)
        }
    )$value
}

# filter junction clusters by browser coordinates
svx_filterJunctionsByRange <- function(jxns, coord, rangeType, chromOnly = FALSE){
    startSpinner(session, message = "filtering junctions to window")
    isWholeGenome <- coord$chromosome == "all"
    if(!isWholeGenome && chromOnly) jxns <- jxns[
        cChrom1 == coord$chromosome & 
        cChrom2 == coord$chromosome
    ]
    switch(
        rangeType,
        center = {
            jxns <- jxns[data.table::between(
                as.numeric(if(isWholeGenome) nodeCenter else refPosCenter), 
                as.numeric(coord$start), 
                as.numeric(coord$end)
            )]
            jxns[ , center := if(isWholeGenome) nodeCenter else refPosCenter]
            jxns
        },
        endpoint = {
            jxns[ , ":="(
                pos1 = if(isWholeGenome) abs(node1) else cRefPos1,
                pos2 = if(isWholeGenome) abs(node2) else cRefPos2
            )]
            jxns[ , ":="(
                pos1In = (isWholeGenome | cChrom1 == coord$chrom) & data.table::between(as.numeric(pos1), as.numeric(coord$start), as.numeric(coord$end)),
                pos2In = (isWholeGenome | cChrom2 == coord$chrom) & data.table::between(as.numeric(pos2), as.numeric(coord$start), as.numeric(coord$end))
            )]
            jxns[pos1In | pos2In]
        },
        jxns
    )   
}

# get all filtered junctions, with (for plots) or without (for nav table) filtering to the coordinate range
svx_getTrackJunctions <- function(track, selectedSources, loadFn, 
                                  coord = NULL, rangeType = NULL, chromOnly = FALSE){
    jxns <- do.call(rbind, c(lapply(names(selectedSources), function(sourceId){
        jxns <- svx_filterJunctionsBySettings(track, sourceId, selectedSources[[sourceId]]$Sample_ID, loadFn)
        cbind(
            jxns,
            sourceId = if(nrow(jxns) == 0) character() else sourceId
        )
    }), list(fill = TRUE))) # tables from different sources have non-identical columns, so must fill
    if(is.null(coord)) jxns else svx_filterJunctionsByRange(jxns, coord, rangeType, chromOnly)
}
