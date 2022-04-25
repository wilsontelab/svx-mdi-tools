#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# objects instantiated here are available to all appStep modules in a session

# get the data.table of called SVs matching the current filter set
getWorkingSvs <- function(settings, sampleSelector){
    assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
    req(assignments)
    req(nrow(assignments) > 0)
    startSpinner(session, 'getWorkingSvs')
    filters <- settings$SV_Filters()
    targetClasses <- SVX$targetClasses[[filters$Target_Class$value]]
    samples <- sampleSelector$selectedSamples()
    isAllSamples <- all(sampleSelector$allSamples() %in% samples)
    setkey(SVX$jxnTypes, name)
    x <- do.call(rbind, lapply(assignments[, unique(Source_ID)], function(sourceId){
        project <- assignments[Source_ID == sourceId][1, Project]
        
        # load table into cache
        svFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "structuralVariants") 

        # apply SV filters
        x <- persistentCache[[svFile]]$data[
            TARGET_CLASS %in% targetClasses &
            JXN_TYPE %in% SVX$jxnTypes[filters$SV_Type$value, code] &
            (JXN_TYPE == "T" | 
             (SV_SIZE >= filters$Min_SV_Size$value & 
              SV_SIZE <= filters$Max_SV_Size$value)) &
            STRAND_COUNT >= filters$Min_Read_Count$value &
            N_SAMPLES >= filters$Min_Samples_With_SV$value &
            N_SAMPLES <= filters$Max_Samples_With_SV$value &
            N_TOTAL >= filters$Min_Source_Molecules$value &
            N_TOTAL <= filters$Max_Source_Molecules$value & 
            N_SPLITS >= filters$Min_Split_Reads$value
        ]

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
    stopSpinner(session, 'getWorkingSvs')
    x[sample.int(.N)] # randomize the list for optimized plotting
}

# get the data.table of molecules that support the currently selected SV
getWorkingMols <- function(sv, sampleSelector){
    req(sv)
    assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
    req(assignments)
    req(nrow(assignments) > 0)
    samples <- unlist(sv$PROJECT_SAMPLES) 
    startSpinner(session, 'getWorkingMols')
    x <- do.call(rbind, lapply(assignments[uniqueId %in% samples, unique(Source_ID)], function(sourceId){
        project <- assignments[Source_ID == sourceId][1, Project]
        molFile <- loadPersistentFile(sourceId = sourceId, contentFileType = "junctionMolecules")
        persistentCache[[molFile]]$data[SV_ID == sv$SV_ID]
    }))
    stopSpinner(session, 'getWorkingSvs')
    x[sample.int(.N)] # randomize the list for optimized plotting
}

# appServer function called after all modules are instantiated
appServer <- function(){


}
