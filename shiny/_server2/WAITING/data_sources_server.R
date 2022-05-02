
#----------------------------------------------------------------------
# ../data_sources_server.R defines data paths and
# loads project and sample lists into server global environment
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# file paths
#----------------------------------------------------------------------
projectsDir  <- paste(serverEnv$DATA_PATH, 'projects', sep="/")
resourcesDir <- paste(serverEnv$DATA_PATH, 'resources', sep="/")
getProjectDir <- function(project){
    paste(projectsDir, project, sep="/")
}
getCaptureTargetsBed <- function(projectInfo){
    file <- paste(projectInfo$CAPTURE_BED, 'bed', sep=".")
    paste(resourcesDir, 'capture_targets', file, sep="/")
}
getExpectedCnvsBed <- function(projectInfo){
    file <- paste(projectInfo$EXPECTED_CNVS_BED, 'bed', sep=".")
    paste(resourcesDir, 'expected_cnvs', file, sep="/")
}
getSamplePrefix <- function(projectInfo, sample){
    paste(getProjectDir(projectInfo$PROJECT), sample, sample, sep="/")
}
getProjectCompareFile <- function(projectInfo, ext='gz'){
    comparePrefix <- paste(getProjectDir(projectInfo$PROJECT), projectInfo$PROJECT, sep="/")
    paste(comparePrefix, projectInfo$GENOME, 'compare', 'structural_variants', ext, sep=".")
}
getSampleGenomePrefix <- function(projectInfo, sample){
    paste(getSamplePrefix(projectInfo, sample), projectInfo$GENOME, sep=".")
}
getExtractFile <- function(projectInfo, sample, fileTypeWithExt){
    paste(getSampleGenomePrefix(projectInfo, sample), 'extract', fileTypeWithExt, sep=".")
}
getSamplePipe <- function(projectInfo, sample, fileType, ext='gz'){
    glob <- paste(getSampleGenomePrefix(projectInfo, sample), 'find', fileType, '*', ext, sep=".")
    cat <- if(ext=='gz') 'zcat' else 'cat'
    pipe(paste(cat, glob))
}
getIndexedFile <- function(projectInfo, sample, fileType, thread){
    paste(getSampleGenomePrefix(projectInfo, sample), 'find', fileType, thread, 'txt', sep=".")
}
getFindRDataFile <- function(projectInfo, sample){
    paste(getSampleGenomePrefix(projectInfo, sample), 'find', 'structural_variants', 'RData', sep=".")  
}
getFindAllNodesFile <- function(projectInfo, sample){
    paste(getSampleGenomePrefix(projectInfo, sample), 'find', 'all_nodes', 'txt', sep=".")  
}
getSampleRDataFile <- function(projectInfo, sample){
    paste(getSampleGenomePrefix(projectInfo, sample), '_server2', serverEnv$SERVER_SUBTYPE, 'RData', sep=".")
}
getSVJunctionImage <- function(project, sample, svId, type="svReads"){
    projectDir <- getProjectDir(project)
    svReadsDir <- file.path(projectDir, sample, 'plots', 'svReads')
    if(!dir.exists(svReadsDir)) dir.create(svReadsDir)
    filename <- paste(sample, type, svId, 'png', sep=".")
    file.path(svReadsDir, filename)
}
