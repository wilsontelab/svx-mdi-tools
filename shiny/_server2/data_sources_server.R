
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

#----------------------------------------------------------------------
# initialize the data sources for the overall app, e.g. all samples
# stays in memory, shared by all user sessions
#----------------------------------------------------------------------
# all available projects that have server ALLOW flag set in project info
# project without this flag are deemed inappropriate for this server interface
alwaysAllow <- list(Distribution=1)
setProjects <- function(){
    reportProgress('setProjects')
    allow <- toupper(paste('allow', serverEnv$SERVER_SUBTYPE, sep="_")) 
    projects <- c('-')
    for(project in list.dirs(projectsDir, full.names=FALSE, recursive=FALSE)){
        infoFile <- paste(getProjectDir(project), 'project.info.txt', sep="/")
        if(file.exists(infoFile)){
            pi <- read.table(infoFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
            allowProject <- pi[pi$name==allow,'value']
            if(length(allowProject)==0) allowProject <- 'FALSE'
            if(!is.null(alwaysAllow[[serverEnv$SERVER_SUBTYPE]])) allowProject <- 'TRUE'
            if(allowProject == 'TRUE') projects <- c(projects, project)             
        }
    }
    assign('projects', projects, envir=.GlobalEnv)
}
# all available samples, stratified by project
.ignoreSamples <- c('plots','svCapture_logs')
setSamples <- function(){
    reportProgress('setSamples')
    samples <- lapply(projects, function(project){
        smps <- list.dirs(getProjectDir(project), full.names=FALSE, recursive=FALSE)
        c('-', smps[smps %notin% .ignoreSamples])
    })
    names(samples) <- projects
    assign('samples', samples, envir=.GlobalEnv)
}

