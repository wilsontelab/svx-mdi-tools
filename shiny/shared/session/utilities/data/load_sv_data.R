#----------------------------------------------------------------------
# functions to load SVs and associated evidence molecules
#----------------------------------------------------------------------

# get the data.table of called SVs matching the current filter set
getFilteredSvs <- function(settings, sampleSelector, isCapture = FALSE){
    assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
    req(assignments)
    req(nrow(assignments) > 0)
    startSpinner(session, 'getFilteredSvs')
    svFilters <- settings$SV_Filters()
    samples <- sampleSelector$selectedSamples()
    isAllSamples <- all(sampleSelector$allSamples() %in% samples)
    setkey(SVX$jxnTypes, name)
    x <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
        project <- assignments[Source_ID == sourceId][1, Project]
        
        # load table into cache
        svFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "structuralVariants") 

        # apply SV filters
        x <- persistentCache[[svFile]]$data[
            JXN_TYPE %in% SVX$jxnTypes[svFilters$SV_Type$value, code] &
            (JXN_TYPE == "T" | 
             (SV_SIZE >= svFilters$Min_SV_Size$value & 
              SV_SIZE <= svFilters$Max_SV_Size$value)) &
            N_SAMPLES >= svFilters$Min_Samples_With_SV$value &
            N_SAMPLES <= svFilters$Max_Samples_With_SV$value &
            N_TOTAL >= svFilters$Min_Source_Molecules$value &
            N_TOTAL <= svFilters$Max_Source_Molecules$value & 
            N_SPLITS >= svFilters$Min_Split_Reads$value
        ]

        # apply capture filters, if applicable
        if(isCapture){
            captureFilters <- if(isCapture) settings$Capture_Filters() else list()
            targetClasses <- SVX$targetClasses[[captureFilters$Target_Class$value]]
            x <- x[
                TARGET_CLASS %in% targetClasses &
                STRAND_COUNT >= captureFilters$Min_Read_Count$value 
            ]  
        }

        # apply sample filters, unless showing all samples
        if(!isAllSamples){
            x[, matchesSamples := FALSE]
            for(sample in samples){
                d <- strsplit(sample, ":")[[1]]
                if(d[1] == project) x[, 
                    matchesSamples := matchesSamples | x[[d[2]]] > 0
                ]
            }
        } else {
            x[, matchesSamples := TRUE]
        }   
        x[, PROJECT_SAMPLES := lapply(strsplit(SAMPLES, ","), function(x) paste(project, x, sep = ":"))]  
        x[matchesSamples == TRUE, .SD, .SDcols = c(names(SVX$find$structural_variants), "PROJECT_SAMPLES") ]
    }))
    stopSpinner(session, 'getFilteredSvs')
    x[sample.int(.N)] # randomize the list for optimized plotting
}

# get a data.table of molecules that support the currently selected SV
#   one SV (one line in the SVs table) may be shared by multiple samples in a source
#   that SV might be listed twice if the selected samples come from multiple sources / pipeline runs
getSVMolecules <- function(sampleSelector, selectedSv){
    assignments <- sampleSelector$selectedAssignments() 
    req(assignments)
    req(nrow(assignments) > 0)
    sv <- selectedSv()
    req(sv)
    startSpinner(session, 'getSVMolecules')    
    samples <- unlist(sv$PROJECT_SAMPLES) 
    sourceId <- assignments[uniqueId %in% samples, unique(Source_ID)] # one SV comes from exactly one source
    molFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "junctionMolecules")
    metadataFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "metadata")
    maxTLen <- max(as.integer(strsplit(as.character(persistentCache[[metadataFile]]$data$MAX_TLENS), "\\s+")[[1]]))
    mols <- persistentCache[[molFile]]$data[SV_ID == sv$SV_ID][sample.int(.N)] # randomize for optimized plotting
    mols[, c('chrom1', 'side1', 'pos1') := unpackNodeNames(NODE_1)]
    mols[, c('chrom2', 'side2', 'pos2') := unpackNodeNames(NODE_2)]
    mols[, isOuterClip := NODE_CLASS == SVX$nodeClasses$OUTER_CLIP]
    stopSpinner(session, 'getSVMolecules')    
    list(
        assignments = assignments,
        sv = sv,
        mols = mols,      
        maxTLen = maxTLen
    )
}

# get the whole genome coverage map

# TODO: need to do this by sample - it is a sample level map, not a project level map

getGenomeCoverage <- function(sampleSelector){
    assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
    req(assignments)
    req(nrow(assignments) > 0)
    startSpinner(session, 'getGenomeCoverage')    
    x <- mapply(function(sourceId, sampleId){
        coverageIndexFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "coverageIndex")
        metadataFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "metadata")

        chromNames <- strsplit(persistentCache[[metadataFile]]$data$CHROMS, "\\s+")[[1]]
        chroms <- seq_along(chromNames)
        names(chroms) <- chromNames

        list(
            chroms = chroms,
            depth = persistentCache[[coverageIndexFile]]$data
        )

    }, assignments$Source_ID, assignments$Sample_ID)
#             Classes 'data.table' and 'data.frame':  46633 obs. of  4 variables:
#  $ chromIndex: int  1 1 1 1 1 1 1 1 1 1 ...
#  $ chunkIndex: int  0 1 2 3 4 5 6 7 8 9 ...
#  $ cumNBreaks: int  0 0 0 0 0 0 0 0 0 0 ...
#  $ coverage  : num  0 0 0 0 0 0 0 0 0 0 ...
    stopSpinner(session, 'getGenomeCoverage')
    x
}
