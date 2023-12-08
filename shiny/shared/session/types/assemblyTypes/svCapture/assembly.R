# ----------------------------------------------------------------------
# top-level support for svCapture assembly loading and parsing
# ----------------------------------------------------------------------
svCapture_getDenominatorColumn <- function(settings){
    x <- "kbp_on_target"
    if(isTruthy(settings$get("SV_Capture","Enforce_Min_Family_Size"))) x <- paste(x, "filtered", sep = "_")
    if(isTruthy(settings$get("SV_Capture","Adjust_Read_Lengths")))     x <- paste(x, "effective", sep = "_")
    x
}
svCapture_loadAssembly <- function(assemblyOptions, settings, rdsFile, ...){
    assemblyCache$get(
        'assembly', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            settings,            
            rdsFile,
            assemblyOptions$internalUseSampleColumns
        ), 
        createFn = function(...) {
            assembly <- readRDS(rdsFile)
            assembly$svs[, edgeType := svx_jxnType_altCodeToX(JXN_TYPE, "code")]
            denominatorColumn <- svCapture_getDenominatorColumn(settings)
            for(col in names(assembly$samples)){
                if(col %in% assemblyOptions$internalUseSampleColumns) next
                values <- assembly$samples[[col]] # ensure that all empty/zero-dose cells have the value "-"
                assembly$samples[[col]][is.na(values) | is.null(values) | values == "" | values == "0"] <- "-"

                ######################## TODO: fix this in samples list and rerun assembly, then delete
                assembly$samples[[col]][values == "bulk"] <- "-"

                values <- assembly$samples[[col]] # remove columns that never vary over all samples
                if(length(unique(values)) == 1) assembly$samples[[col]] <- NULL
            }
            assembly$samples$coverage <- as.integer(assembly$samples[[denominatorColumn]] / assembly$samples$sumTargetLens * 1000)
            assembly
        }
    )$value
}
svCapture_getGroups <- function(assemblyOptions, groupedProjectSamples, groupingCols, input, ...){
    assemblyCache$get(
        'groups', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            groupedProjectSamples,
            groupingCols,
            input$Data_Types
        ), 
        createFn = function(...) {
            gps <- groupedProjectSamples
            gps$nSvs <- 0
            gps$svFrequency <- 0
            for(svType in input$Data_Types) {
                nSvs <- as.integer(gps[[svType]])
                gps$nSvs <- gps$nSvs + nSvs
                gps$svFrequency <- gps$svFrequency + nSvs / gps$coverage
            }
            groups <- gps[, .(
                nProjects = length(unique(project)), # must set nProject, nSamples
                nSamples = .N,
                coverage = sum(coverage),
                nSvs = sum(nSvs),
                meanFrequency = round(mean(svFrequency), 4),
                sdFrequency = round(sd(svFrequency), 4),
                svFrequencies = list(svFrequency),
                projects = list(project)
            ), by = groupingCols]
            setAssemblyGroupLabels(groups, groupingCols)
        }
    )$value
}
#----------------------------------------------------------------------
# SVs found in the working samples, to populate distribution plots
#----------------------------------------------------------------------
matchingSvs <- function(assemblyOptions, isProcessingData, assembly, groupedProjectSamples, groupingCols, input){
    req(isProcessingData())
    assembly <- assembly()
    gps <- groupedProjectSamples() 
    groupingCols <- groupingCols()   
    req(assembly, gps, groupingCols, input$Data_Types)
    assemblyCache$get(
        'matchingSvs', 
        permanent = TRUE,
        from = "ram",
        create = assemblyOptions$cacheCreateLevel,
        keyObject = list(
            assembly$svs,
            gps,
            groupingCols,
            input$Data_Types
        ), 
        createFn = function(...) {
            startSpinner(session, message = "collecting matching SVs")
            x <- gps[, {
                assembly$svs[
                    PROJECT == project & 
                    grepl(paste0(",", sample, ","), SAMPLES) & 
                    edgeType %in% svx_jxnType_longNameToX(input$Data_Types, "code")
                ]
            }, by = c("project", "sample", groupingCols)]
            x <- setAssemblyGroupLabels(x, groupingCols)
            stopSpinner(session)
            x
        }
    )$value
}
