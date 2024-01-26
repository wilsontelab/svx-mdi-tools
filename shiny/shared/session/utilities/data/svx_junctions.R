#----------------------------------------------------------------------
# handle unique junction loading and filtering
# these are generic functions that apply to potentially any SVX app
# expects:
#   sessionCache
#   loadFn(targetId) that returns a data.table with appropriately named display and filter columns
#----------------------------------------------------------------------
svxLoadCreate <- "asNeeded" # for convenience when debugging load functions

# load all unique junctions for a specific targetId (not filtered by region, type, or sample yet)
svx_loadJunctions <- function(targetId, loadFn){
    req(targetId)
    sessionCache$get(
        'svx_loadJunctions', 
        key = targetId, 
        permanent = TRUE, 
        from = "ram", 
        create = svxLoadCreate, 
        createFn = function(...) {
            startSpinner(session, message = "loading junctions")
            jxns <- loadFn(targetId)
            getCenter <- function(p1, p2) pmin(p1, p2) + abs(p2 - p1) / 2
            jxns[, samples :=  paste0(",", samples, ",")] # for simpler matching of single samples
            jxns[, ":="(
                nodeCenter   = getCenter(abs(node1), abs(node2)),
                refPosCenter = getCenter(cRefPos1,   cRefPos2),
                size = abs(abs(node2) - abs(node1)) # eventSize has the chrom-level size, where translocation size == 0
            )]   
            tryCatch({ # if no chromosomesFile available (e.g., svAmplicon), loadFn must set cChrom1/2
                chroms <- readRDS(getSourceFilePath(targetId, "chromosomesFile")) 
                jxns[, ":="(
                    cChrom1 = unlist(chroms$revChromIndex[cChromIndex1]),
                    cChrom2 = unlist(chroms$revChromIndex[cChromIndex2])
                )]                    
            }, error = function(e) {
                print(e)
                NULL
            })
            jxns
        }
    )$value
}

# filter unique junctions by sample, e.g., specified by browser items list
# this applies when target items are samples, but junctions are loaded from a parent sourceId carrying multiple samples
svx_filterJunctionsBySample <- function(targetId, samples_, loadFn){
    sessionCache$get(
        'svx_filterJunctionsBySample', 
        keyObject = list(targetId = targetId, samples = samples_), 
        permanent = TRUE, 
        from = "ram",
        create = svxLoadCreate, 
        createFn = function(...) {
            jxns <- svx_loadJunctions(targetId, loadFn)
            if(is.null(samples_)) return(jxns) # e.g., when track items are amplicons
            startSpinner(session, message = "filtering junctions")
            I <- FALSE
            samples_ <- if(startsWith(samples_[1], ",")) samples_ else paste0(",", samples_, ",")
            for(sample_ in samples_) I <- I | jxns[, grepl(sample_, samples)]
            jxns[I]
        }
    )$value
}

# when in use in app, add a CN column to CNV junctions
svx_setCnvJxnCN <- function(jxns, targetId, samples_, cnData){
    if(is.na(cnData$key)) return(jxns)
    sessionCache$get(
        'svx_setCnvJxnCN', 
        keyObject = list(targetId = targetId, samples = samples_, cnDataKey = cnData$key), 
        permanent = FALSE, 
        from = "ram",
        create = svxLoadCreate, 
        createFn = function(...) {
            startSpinner(session, message = "adding junction CN")
            merge(
                jxns,
                cnData$value[, .SD, .SDcols = c("SV_ID", "outerFlankCN", "innerFlankCN", "flankCNC", "junctionCN")],
                by = "SV_ID",
                all.x = TRUE
            )
        }
    )$value
}

# filter junction clusters based on track or other settings
svx_applyJunctionFilters <- function(jxns, settings){

    startSpinner(session, message = paste("filtering junctions"))

    filters <- if(is.null(settings$Filters)) list() else settings$Filters()
    filters <- lapply(names(svx_filterDefaults), function(filter){
        if(is.null(filters[[filter]])) svx_filterDefaults[[filter]]
        else if(!is.null(filters[[filter]]$selected)) filters[[filter]]$selected else filters[[filter]]$value
    })
    names(filters) <- names(svx_filterDefaults)

    if(filters$Min_SV_Size > 1) jxns <- jxns[size >= filters$Min_SV_Size]
    if(filters$Max_SV_Size > 0) jxns <- jxns[size <= filters$Max_SV_Size] #  & edgeType != svx_edgeTypes$TRANSLOCATION

    jxns <- jxns[insertSize >= filters$Min_Insert_Size] # can be negative
    jxns <- jxns[insertSize <= filters$Max_Insert_Size]

    if(filters$Min_Source_Molecules > 1) jxns <- jxns[nInstances >= filters$Min_Source_Molecules]
    if(filters$Max_Source_Molecules > 0) jxns <- jxns[nInstances <= filters$Max_Source_Molecules]

    if(filters$Min_Sequenced_Molecules > 0) jxns <- jxns[nSequenced >= filters$Min_Sequenced_Molecules]
    if(filters$Max_Linked_Junctions > 0) jxns <- jxns[nLinkedJunctions <= filters$Max_Linked_Junctions]

    if(filters$Min_Samples_With_SV > 1) jxns <- jxns[nSamples >= filters$Min_Samples_With_SV]
    if(filters$Max_Samples_With_SV > 0) jxns <- jxns[nSamples <= filters$Max_Samples_With_SV]
    if(filters$Unique_To_Sample != "show all SVs") jxns <- jxns[samples == paste0(",", filters$Unique_To_Sample, ",")]

    if(filters$Show_ChrM != "always"){
        hasChrM <- jxns[, cChrom1 == "chrM" | cChrom2 == "chrM"]
        if(filters$Show_ChrM == "translocations_only") jxns <- jxns[hasChrM == FALSE | xor(cChrom1 == "chrM", cChrom2 == "chrM")]
        else if(filters$Show_ChrM == "never")          jxns <- jxns[hasChrM == FALSE]
    }

    if(filters$Min_Map_Quality > 0)jxns <- jxns[mapQ >= filters$Min_Map_Quality] 
    if(filters$Min_Flank_Length > 0) jxns <- jxns[flankLength >= filters$Min_Flank_Length] 
    if(filters$Min_Flank_CNC > 0 && "flankCNC" %in% names(jxns)) jxns <- jxns[!is.na(flankCNC) & abs(flankCNC) >= filters$Min_Flank_CNC] 

    if(length(filters$SV_Type) > 0) {
        edgeTypes <- svx_jxnType_nameToX(filters$SV_Type, "code")
        jxns <- jxns[edgeType %in% edgeTypes]
    }

    jxns
}
svx_filterJunctionsBySettings <- function(track, targetId, samples, loadFn, family, cnData = NULL){
    jxns <- sessionCache$get(
        'svx_filterJunctionsBySettings', 
        keyObject = list(targetId = targetId, samples = samples, cnDataKey = cnData$key, settings = track$settings$all()), 
        permanent = FALSE,
        from = "ram", 
        create = svxLoadCreate, 
        createFn = function(...) {
            jxns <- svx_filterJunctionsBySample(targetId, samples, loadFn) %>% 
                    svx_setCnvJxnCN(targetId, samples, cnData) %>%
                    svx_applyJunctionFilters(track$settings)
            stopSpinner(session)
            jxns %>% svx_setJunctionPointColors(track, family) %>% svx_setJunctionPointSizes(track)
        }
    )$value
}

# filter junction clusters by browser coordinates
svx_filterJunctionsByRange <- function(jxns, coord, rangeType, chromOnly = FALSE){
    startSpinner(session, message = "filtering junctions to window")
    isWholeGenome <- coord$chromosome == "all"
    isProperChrom <- isProperChromosome(coord$chromosome)
    if(isProperChrom && chromOnly) jxns <- jxns[
        cChrom1 == coord$chromosome & 
        cChrom2 == coord$chromosome
    ]
    switch(
        rangeType,
        centerAll = {
            jxns[ , center := if(isProperChrom) refPosCenter else nodeCenter]
            jxns
        },
        center = {
            jxns <- jxns[data.table::between(
                as.numeric(if(isProperChrom) refPosCenter else nodeCenter), 
                as.numeric(coord$start), 
                as.numeric(coord$end)
            )]
            jxns[ , center := if(isProperChrom) refPosCenter else nodeCenter]
            jxns
        },
        endpoint = {
            jxns[ , ":="(
                pos1 = if(isProperChrom) cRefPos1 else abs(node1),
                pos2 = if(isProperChrom) cRefPos2 else abs(node2)
            )]
            jxns[ , ":="(
                pos1In = (!isProperChrom | cChrom1 == coord$chrom) & data.table::between(as.numeric(pos1), as.numeric(coord$start), as.numeric(coord$end)),
                pos2In = (!isProperChrom | cChrom2 == coord$chrom) & data.table::between(as.numeric(pos2), as.numeric(coord$start), as.numeric(coord$end))
            )]
            jxns[pos1In | pos2In]
        },
        jxns
    )   
}

# get all filtered junctions, with (for plots) or without (for nav table) filtering to the coordinate range
#   when isMultiSample == TRUE  (e.g., svPore), selected targets are samples, junctions are loaded from the parent sources then filtered by sample
#   when isMultiSample == FALSE (e.g., svAmplicon), selected targets are end-level (e.g., amplicon keys), and junctions are loaded directly from there
svx_getTrackJunctions <- function(track, selectedTargets, loadFn, 
                                  coord = NULL, rangeType = NULL, chromOnly = FALSE,
                                  isMultiSample = TRUE, family = "Points"){ 
    jxns <- do.call(rbind, c(lapply(names(selectedTargets), function(targetId){
        jxns <- svx_filterJunctionsBySettings(
            track, 
            targetId, 
            if(isMultiSample) selectedTargets[[targetId]]$Sample_ID else NULL,
            loadFn,
            family,
            cnData = svx_getCnvJxnNormalizedCN(isMultiSample, targetId)
        )
        cbind(
            jxns,
            targetId = if(nrow(jxns) == 0) character() else targetId
        )
    }), list(fill = TRUE))) # tables from different sources have non-identical columns, so must fill
    if(is.null(coord)) jxns else svx_filterJunctionsByRange(jxns, coord, rangeType, chromOnly)
}

# support composite genome colors as purple inter-genomic translocation dots/lines
svx_colorCompositeGenomes <- function(x, reference){
    if(objectHasData(listCompositeGenomes(reference))){
        genome1 <- sapply(strsplit(x$cChrom1, "_"), function(x) x[2])
        genome2 <- sapply(strsplit(x$cChrom2, "_"), function(x) x[2])
        x[genome1 != genome2, ":="(
            color = CONSTANTS$plotlyColors$purple,
            edgeType = svx_edgeTypes$INTERGENOME,
            interGenome = TRUE
        )]
    }
    x
}
