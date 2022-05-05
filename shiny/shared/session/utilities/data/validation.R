#----------------------------------------------------------------------
# data validation functions
#----------------------------------------------------------------------

# reject a sampleSet if samples present in manifest but not in processed data file
validateSampleAssignments <- function(sampleSet){
    checkSampleDataExists(sampleSet, function(sourceId){
        path <- getSourceFilePath(sourceId, 'sampleSummary')
        sampleSummary <- fread(path)
        knownSampleIds <- c(sampleSummary$sampleName, sampleSummary$libraryName)        
    })
}
